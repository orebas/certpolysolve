// Reusable interpreter slots: exact operation order, bounded live storage,
// early outputs/dead nodes, and cleanup plus reuse after every operation throws.
#include "numeric/evaluate/program_arithmetic.hpp"
#include "numeric/evaluate/program_interpreter.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <algorithm>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>

using namespace certpolysolve::detail;
namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
static_assert(!std::is_constructible_v<PreparedProgram, Program &&>);
static_assert(!std::is_constructible_v<PreparedProgram, const Program &&>);

struct Trace {
    size_t live = 0, peak = 0, calls = 0;
    std::optional<size_t> fail;
    std::vector<ProgramOp> operations;
};
struct TrackedArithmetic {
    struct Value {
        DualBallArithmetic::Value balls;
        Trace *trace;
        Value(DualBallArithmetic::Value value, Trace &t) : balls(std::move(value)), trace(&t) {
            t.peak = std::max(t.peak, ++t.live);
        }
        Value(Value &&other) noexcept : balls(std::move(other.balls)), trace(std::exchange(other.trace, nullptr)) {}
        Value(const Value &) = delete;
        ~Value() { if (trace) --trace->live; }
    };
    DualBallArithmetic arithmetic;
    Trace &trace;
    TrackedArithmetic(const acb_t z, slong bits, Trace &t) : arithmetic(z, bits), trace(t) {}
    template <class Call> Value apply(ProgramOp op, Call call) {
        const size_t position = trace.calls++;
        if (trace.fail == position) throw std::bad_alloc();
        trace.operations.push_back(op);
        return Value(call(), trace);
    }
    Value variable() { return apply(ProgramOp::Variable, [&] { return arithmetic.variable(); }); }
    Value constant(const Coefficient &c) { return apply(ProgramOp::Constant, [&] { return arithmetic.constant(c); }); }
    Value add(const Value &a, const Value &b) {
        return apply(ProgramOp::Add, [&] { return arithmetic.add(a.balls, b.balls); });
    }
    Value subtract(const Value &a, const Value &b) {
        return apply(ProgramOp::Subtract, [&] { return arithmetic.subtract(a.balls, b.balls); });
    }
    Value multiply(const Value &a, const Value &b) {
        return apply(ProgramOp::Multiply, [&] { return arithmetic.multiply(a.balls, b.balls); });
    }
    Value negate(const Value &a) { return apply(ProgramOp::Negate, [&] { return arithmetic.negate(a.balls); }); }
    Value power(const Value &a, uint64_t exponent) {
        return apply(ProgramOp::Power, [&] { return arithmetic.power(a.balls, exponent); });
    }
};

void exercise(const Program &program, size_t max_slots, bool inject) {
    const PreparedProgram prepared(program);
    require(prepared.slots.count <= max_slots, "workspace grows with dead or expired nodes");
    std::vector<std::optional<TrackedArithmetic::Value>> values;
    std::vector<ProgramOp> order;
    for (const auto &node : program.nodes) order.push_back(node.op);
    Trace trace;
    const auto empty = [&] {
        return std::all_of(values.begin(), values.end(), [](const auto &v) { return !v; });
    };
    for (slong bits : {slong{64}, slong{512}, slong{64}}) {
        ComplexBall z;
        acb_set_d_d(z.value(), 0.375, -0.3125);
        arb_add_error_2exp_si(acb_realref(z.value()), -50);
        DualBallArithmetic reference(z.value(), bits);
        const auto expected = run_program(program, reference);
        const auto check = [&] {
            trace.calls = trace.peak = 0;
            trace.operations.clear();
            {
                TrackedArithmetic arithmetic(z.value(), bits, trace);
                const auto got = run_prepared(prepared, arithmetic, values);
                require(acb_equal(got.balls.value.value(), expected.value.value()) &&
                        acb_equal(got.balls.derivative.value(), expected.derivative.value()),
                        "workspace changes the value or derivative balls");
                require(trace.operations == order, "workspace changes arithmetic order or omits dead nodes");
                require(trace.live == 1 && empty(), "temporaries survive a completed evaluation");
            }
            require(trace.live == 0 && trace.peak <= max_slots, "live value ownership or peak is wrong");
            require(values.size() <= max_slots, "fixed storage follows node count instead of liveness");
        };
        check();
        if (!inject) continue;
        for (size_t failure = 0; failure < program.nodes.size(); ++failure) {
            trace.calls = trace.peak = 0;
            trace.operations.clear();
            trace.fail = failure;
            bool threw = false;
            try {
                TrackedArithmetic arithmetic(z.value(), bits, trace);
                (void)run_prepared(prepared, arithmetic, values);
            } catch (const std::bad_alloc &) { threw = true; }
            require(threw && trace.calls == failure + 1, "injection missed its arithmetic operation");
            require(trace.live == 0 && empty(), "temporaries survive an interrupted evaluation");
            trace.fail.reset();
            check();
        }
    }
}

Program recurrence() {
    Program program = parse_program("PSPROGRAM 1\nx\nc 1\nc 2\nmul 2 0\nout 0\n");
    uint64_t previous = 1, current = 0;
    for (size_t i = 1; i < 512; ++i) {
        const auto product = program.append({ProgramOp::Multiply, 3, current, 0});
        const auto next = program.append({ProgramOp::Subtract, product, previous, 0});
        previous = current;
        current = next;
    }
    program.output = current;
    return program;
}
} // namespace

int main() {
    certpolysolve::ThreadScope cleanup;
    try {
        // Output is read by a later unused node; repeated operands must release once.
        exercise(parse_program("PSPROGRAM 1\nx\nc 3\nmul 0 0\nadd 2 1\nneg 3\npow 1 2\nout 3\n"), 4, true);
        // A constant/variable with no readers must release its slot immediately.
        exercise(parse_program("PSPROGRAM 1\nx\nc 7\nx\nx\nc 5\nout 0\n"), 2, true);
        exercise(parse_program("PSPROGRAM 1\nc 3\nout 0\n"), 1, true);
        exercise(recurrence(), 6, false);
        std::cout << "program workspace order, peak storage, exceptions and reuse pass\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
