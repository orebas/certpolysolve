#include "api_test_support.hpp"
#include "exact_winding_api.hpp"
#include "numeric/evaluate/program_machine.hpp"
#include <bit>
#include <cstdlib>
#include <iostream>
#include <new>
#include <type_traits>

// Deliberately inject failures into ordinary C++ allocations only. FLINT's
// allocator failure policy is independent and is not claimed to be covered.
static thread_local long remaining = -1;
static thread_local long live_allocations = 0;
void *operator new(std::size_t size) {
    if (remaining == 0) { remaining = -1; throw std::bad_alloc(); }
    if (remaining > 0) --remaining;
    if (auto *p = std::malloc(size ? size : 1)) { ++live_allocations; return p; }
    throw std::bad_alloc();
}
void *operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void *p) noexcept { if (p) --live_allocations; std::free(p); }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete(void *p, std::size_t) noexcept { ::operator delete(p); }
void operator delete[](void *p, std::size_t) noexcept { ::operator delete(p); }

// The C++ wrapper must own a successful C result before its next allocation.
// Check every allocation position through completion, including wrapper work
// after the C operation returns, and require no ordinary C++ allocation to leak.
template <class Operation>
void cpp_operation(const char *name, Operation operation, const certpolysolve::Polynomial &expected) {
    unsigned failures = 0, successes = 0;
    for (long fault = 0; fault < 128; ++fault) {
        const long before = live_allocations;
        remaining = fault;
        try {
            const auto result = operation();
            remaining = -1;
            if (!result.ok() || !(result == expected))
                throw std::runtime_error(std::string(name) + ": wrong result");
            ++successes;
        } catch (const std::bad_alloc &) {
            remaining = -1;
            ++failures;
        }
        const bool balanced = live_allocations == before; // before constructing the failure message
        if (!balanced)
            throw std::runtime_error(std::string(name) + ": leaked allocation at fault " + std::to_string(fault));
    }
    if (!failures || !successes)
        throw std::runtime_error(std::string(name) + ": did not cover failure and success");
    std::cout << name << ": " << failures << " allocation failures, " << successes << " completed operations, no leaks\n";
}

void cpp_arithmetic() {
    const certpolysolve::Polynomial a{1, 2, 1}, b{-1, 1};
    const std::vector<certpolysolve::Number> roots{1, -1, 2};
    cpp_operation("add", [&] { return a + b; }, {0, 3, 1});
    cpp_operation("subtract", [&] { return a - b; }, {2, 1, 1});
    cpp_operation("multiply", [&] { return a * b; }, {-1, -1, 1, 1});
    cpp_operation("derivative", [&] { return a.derivative(); }, {2, 2});
    cpp_operation("power", [&] { return pow(b, 3); }, {-1, 3, -3, 1});
    cpp_operation("negate", [&] { return -a; }, {-1, -2, -1});
    cpp_operation("from_roots", [&] { return certpolysolve::Polynomial::from_roots(roots); }, {2, -1, -2, 1});
}

// Buffer growth can fail after only some live slots have grown. The same owner
// must remain usable for both its old extent and a subsequent growth attempt.
// Exercise every allocation position until the first uninterrupted evaluation;
// the test does not depend on the number or layout of the owned buffers.
template <class R> void program_buffers() {
    using namespace certpolysolve::detail;
    using Machine = ProgramMachine<R>;
    const auto program = parse_program("PSPROGRAM 1\nx\nc 1\nmul 0 0\nadd 2 1\nmul 3 2\nsub 4 0\nout 5\n");
    Machine reference(program);
    std::vector<typename Machine::C> points(257);
    std::vector<typename Machine::Evaluation> expected(points.size()), actual(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        points[i] = {R(0.25 + 0.001 * static_cast<double>(i)), R(-0.125)};
        expected[i] = reference.evaluate_one(points[i]);
    }
    const auto equal = [](R a, R b) {
        const auto bits = [](double x) { return std::bit_cast<uint64_t>(x); };
        if constexpr (std::is_same_v<R, double>) return bits(a) == bits(b);
        else return bits(a.m) == bits(b.m) && a.e == b.e;
    };
    const auto check = [&](const auto &a, const auto &b) {
        if (a.finite != b.finite || !equal(a.value.re, b.value.re) || !equal(a.value.im, b.value.im) ||
            !equal(a.derivative.re, b.derivative.re) || !equal(a.derivative.im, b.derivative.im) ||
            !equal(a.noise, b.noise) || !equal(a.noise_derivative, b.noise_derivative))
            throw std::runtime_error("program buffer allocation failure corrupted a later evaluation");
    };
    for (bool warmed : {false, true}) {
        unsigned failures = 0;
        bool completed = false;
        for (long fault = 0; fault < 128 && !completed; ++fault) {
            const long before = live_allocations;
            {
                Machine reused(program);
                const auto input = std::span<const typename Machine::C>(points);
                const auto output = std::span<typename Machine::Evaluation>(actual);
                if (warmed) reused.evaluate(input.first(3), output.first(3), 3);
                remaining = fault;
                try {
                    reused.evaluate(input, output, points.size());
                    remaining = -1;
                    completed = true;
                } catch (const std::bad_alloc &) {
                    remaining = -1;
                    ++failures;
                }
                reused.evaluate(input.first(3), output.first(3), 3);
                for (size_t i = 0; i < 3; ++i) check(actual[i], expected[i]);
                reused.evaluate(input, output, points.size());
                for (size_t i = 0; i < points.size(); ++i) check(actual[i], expected[i]);
                Machine copied(reused);
                Machine moved(std::move(copied));
                check(moved.evaluate_one(points[0]), expected[0]);
            }
            if (live_allocations != before)
                throw std::runtime_error("program buffer failure/recovery leaked an allocation");
        }
        if (!failures || !completed)
            throw std::runtime_error("program buffer fault injection missed failure or success");
        std::cout << "program buffers (" << sizeof(R) << " bytes, warmed=" << warmed << "): "
                  << failures << " failures recovered, no leaks\n";
    }
}

// A solve interrupted by an allocation failure returns CPS_OUT_OF_MEMORY and publishes no result; every other
// solve's published counts are right. (Until 18 October 2026 this injected the failures into a session's
// advances and checked that the session refused reuse; the 0.5 API has one-shot solves only.)
int main(int argc, char **argv) {
    certpolysolve::ThreadScope lifetime;
    try {
        cpp_arithmetic();
        if (argc == 2 && std::string_view(argv[1]) == "cpp-only") return 0;
        program_buffers<double>();
        program_buffers<certpolysolve::detail::machine::BlockDouble>();
        if (argc == 2 && std::string_view(argv[1]) == "program-only") return 0;
        using namespace reference;
        const std::vector<Complex> known{{-2, 1}, {Rational(1)/3, Rational(1)/5}, {3, -2}};
        const auto polynomial = from_roots(known);
        auto input_poly = input(polynomial);
        unsigned failures = 0, solved = 0;
        for (auto method : {certpolysolve::Method::automatic, certpolysolve::Method::classical, certpolysolve::Method::secular,
                            certpolysolve::Method::piecewise}) {
            // C options: the solve below is the C call itself, so that an injected failure reaches it as
            // CPS_OUT_OF_MEMORY and not as an exception from the header.
            const auto options = certpolysolve_test::raw_options();
            require(cps_options_set_method(options.get(), static_cast<int>(method)) == CPS_OK, "method not set");
            require(cps_options_set_goal(options.get(), CPS_ROOT_APPROXIMATION) == CPS_OK, "goal not set");
            for (long fault = 0; fault < 192; ++fault) {
                cps_result *out = nullptr;
                remaining = fault;
                const auto status = cps_solve(input_poly.get(), options.get(), &out);
                remaining = -1;
                if (status == CPS_OUT_OF_MEMORY) {
                    ++failures;
                    require(out == nullptr, "allocation failure published a result");
                    continue;
                }
                require(status == CPS_OK, "unexpected injected-allocation error");
                ++solved;
                const auto result = certpolysolve_test::from_c(CPS_OK, out);
                require(result.goal_met(), "a solve the allocation failure did not reach failed");
                uint64_t count = 0;
                for (const auto &root : result.roots()) {
                    require(rectangle_count(polynomial, root) == root.count,
                            "allocation failure corrupted a published count");
                    check_accuracy(root, 53);
                    count += root.count;
                }
                require(count == known.size(), "allocation failure lost roots");
            }
        }
        require(failures > 100, "fault injection missed the solves' allocations");
        std::cout << failures << " allocation failures refused; " << solved << " solves past the injected fault\n";
        return 0;
    } catch (const std::exception &error) {
        remaining = -1;
        std::cerr << error.what() << '\n';
        return 1;
    }
}
