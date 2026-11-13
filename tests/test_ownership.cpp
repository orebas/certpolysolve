#include <certpolysolve/certpolysolve.hpp>
#include "numeric/evaluate/evaluator.hpp"
#include <concepts>
#include <iostream>
#include <type_traits>
#include <vector>

namespace d = certpolysolve::detail;

template<class T> concept TemporaryValue = requires(T &&x) { std::move(x).value(); };
template<class T> concept TemporaryData = requires(T &&x) { std::move(x).data(); };
template<class T> concept TemporaryView = requires(T &&x) { std::move(x).view(); };
template<class T> concept TemporaryElement = requires(T &&x) { std::move(x)[0]; };
template<class T> concept TemporaryCoefficients = requires(T &&x) { std::move(x).coefficients(); };
static_assert(!TemporaryValue<d::Integer> && !TemporaryValue<d::Rational> &&
              !TemporaryValue<d::Magnitude> && !TemporaryValue<d::RealBall> &&
              !TemporaryValue<d::Float> && !TemporaryValue<d::ComplexBall> &&
              !TemporaryValue<d::RealPolynomial> && !TemporaryValue<const d::RealPolynomial> &&
              !TemporaryValue<d::ComplexPolynomial> && !TemporaryValue<d::IntegerPolynomial> &&
              !TemporaryValue<d::Factorization> && !TemporaryValue<d::PolynomialSnapshot>);
static_assert(!TemporaryData<d::RootVector> && !TemporaryData<d::RealVector>);
static_assert(!TemporaryView<d::RootVector> && !TemporaryElement<d::RootVector> &&
              !TemporaryElement<d::RealVector>);
static_assert(std::same_as<decltype(std::declval<const d::ComplexBall &>().value()), acb_srcptr>);
static_assert(std::same_as<decltype(std::declval<const d::RootVector &>().data()), acb_srcptr>);
static_assert(std::same_as<decltype(std::declval<const d::RootVector &>()[0]), acb_srcptr>);
static_assert(std::same_as<decltype(std::declval<const d::RealVector &>()[0]), arb_srcptr>);
static_assert(!std::is_copy_constructible_v<d::RootVector> && !std::is_copy_constructible_v<d::RealVector>);
static_assert(!std::is_copy_constructible_v<d::RealPolynomial> && !std::is_copy_assignable_v<d::RealPolynomial>);
static_assert(std::is_nothrow_move_constructible_v<d::RealPolynomial> &&
              std::is_nothrow_move_assignable_v<d::RealPolynomial>);
static_assert(std::same_as<decltype(std::declval<d::RealPolynomial &>().value()), arb_poly_struct *>);
static_assert(std::same_as<decltype(std::declval<const d::RealPolynomial &>().value()), const arb_poly_struct *>);
static_assert(!TemporaryCoefficients<d::ComplexPolynomial> && !TemporaryCoefficients<const d::ComplexPolynomial> &&
              !TemporaryCoefficients<d::PolynomialSnapshot> && !TemporaryCoefficients<const d::PolynomialSnapshot>);
static_assert(std::is_trivially_copyable_v<d::ConstComplexBallView>);
static_assert(std::same_as<decltype(std::declval<const d::ComplexPolynomial &>().coefficients()[0]), acb_srcptr>);
static_assert(std::same_as<decltype(std::declval<const d::PolynomialSnapshot &>().coefficients()[0]), acb_srcptr>);

static void require(bool ok, const char *message) {
    if (!ok) throw std::runtime_error(message);
}
template<class Exception, class F> static void rejects(F &&f, const char *message) {
    try { f(); } catch (const Exception &) { return; }
    throw std::runtime_error(message);
}
static void coefficients(d::ComplexPolynomial &p) {
    acb_poly_set_coeff_si(p.value(), 0, 1);
    acb_poly_set_coeff_si(p.value(), 1, 2);
    acb_poly_set_coeff_si(p.value(), 2, 3);
}
static void check_evaluator(const d::Evaluator &evaluator) {
    d::ComplexBall z, value, derivative, expected_value, expected_derivative;
    acb_set_si(z.value(), 2);
    arb_one(acb_imagref(z.value()));
    acb_set_si(expected_value.value(), 14);
    arb_set_si(acb_imagref(expected_value.value()), 14);
    acb_set_si(expected_derivative.value(), 14);
    arb_set_si(acb_imagref(expected_derivative.value()), 6);
    require(evaluator.degree() == 2, "snapshot degree changed");
    require(evaluator.evaluate(value.value(), derivative.value(), z.value(), 128), "snapshot evaluation failed");
    require(acb_contains(value.value(), expected_value.value()) &&
            acb_contains(derivative.value(), expected_derivative.value()),
            "snapshot no longer encloses 1+2z+3z^2 and its derivative at 2+i");
}

template<class Factory> static void owned_evaluator(Factory factory) {
    std::unique_ptr<d::Evaluator> retained;
    {
        d::ComplexPolynomial original;
        coefficients(original);
        retained = factory(original);
        check_evaluator(*retained);
        acb_poly_set_coeff_si(original.value(), 0, 99);
        check_evaluator(*retained);
        // Changing degree and destroying the caller's storage cannot invalidate
        // an evaluator that retains its own immutable coefficient snapshot.
        acb_poly_zero(original.value());
    }
    check_evaluator(*retained);
}
static void snapshots() {
    owned_evaluator([](const d::ComplexPolynomial &p) {
        return std::make_unique<d::CoefficientEvaluator>(p.value());
    });
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The local Taylor cache (the research build's Series test engine's evaluator).
    owned_evaluator([](const d::ComplexPolynomial &p) {
        return std::make_unique<d::LocalPolynomialEvaluator>(p.value());
    });
#endif
    d::ComplexPolynomial original;
    coefficients(original);
    d::PolynomialSnapshot snapshot(std::move(original));
    original = d::ComplexPolynomial{};
    acb_poly_set_coeff_si(original.value(), 0, 99);
    d::CoefficientEvaluator a(snapshot);
    check_evaluator(a);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    d::LocalPolynomialEvaluator b(snapshot);
    check_evaluator(b);
#endif
    auto retained = std::move(snapshot);
    // NOLINTNEXTLINE(bugprone-use-after-move): deliberately test the checked empty-state contract.
    rejects<std::logic_error>([&] { (void)snapshot.value(); }, "moved snapshot lent a backend pointer");
    d::CoefficientEvaluator still_valid(retained);
    check_evaluator(still_valid);
    rejects<std::invalid_argument>([] { d::PolynomialSnapshot invalid(nullptr); }, "null snapshot accepted");
}
static void real_polynomials() {
    d::RealPolynomial moved;
    {
        d::RealPolynomial original;
        arb_poly_set_coeff_si(original.value(), 2, 7);
        auto *coefficients = original.value()->coeffs;
        d::RealPolynomial constructed(std::move(original));
        require(arb_poly_is_zero(original.value()), "move construction did not empty the source");
        require(constructed.value()->coeffs == coefficients, "move construction copied polynomial storage");
        arb_poly_set_coeff_si(moved.value(), 0, 11);
        moved = std::move(constructed);
        require(moved.value()->coeffs == coefficients, "move assignment copied polynomial storage");
        d::RealBall previous;
        arb_poly_get_coeff_arb(previous.value(), constructed.value(), 0);
        require(arb_equal_si(previous.value(), 11), "move assignment did not retain the old destination in its source");
    }
    d::RealBall coefficient;
    arb_poly_get_coeff_arb(coefficient.value(), moved.value(), 2);
    require(arb_equal_si(coefficient.value(), 7), "destroying the moved source invalidated the polynomial");
    auto &alias = moved;
    moved = std::move(alias);
    arb_poly_get_coeff_arb(coefficient.value(), moved.value(), 2);
    require(arb_equal_si(coefficient.value(), 7), "self-move changed the polynomial");
    // The seed pass holds owners in vectors. Relocation must transfer storage without deep copies.
    std::vector<d::RealPolynomial> blocks;
    blocks.push_back(std::move(moved));
    const auto *coefficients = blocks[0].value()->coeffs;
    blocks.reserve(blocks.capacity() + 1);
    require(blocks[0].value()->coeffs == coefficients, "vector relocation copied polynomial storage");
    arb_poly_get_coeff_arb(coefficient.value(), blocks[0].value(), 2);
    require(arb_equal_si(coefficient.value(), 7), "vector relocation changed the polynomial");
}
static void coefficient_views() {
    d::ComplexPolynomial original;
    require(original.coefficients().size() == 0, "empty polynomial supplied a nonempty coefficient view");
    rejects<std::out_of_range>([&] { (void)original.coefficients()[0]; }, "empty coefficient view accepted indexing");
    coefficients(original);
    const auto view = original.coefficients();
    require(view.size() == 3 && acb_equal_si(view[0], 1) && acb_equal_si(view[1], 2) && acb_equal_si(view[2], 3),
            "coefficient view changed extent, order or values");
    require(view[0] == original.value()->coeffs, "coefficient view copied storage");
    rejects<std::out_of_range>([&] { (void)view[-1]; }, "negative coefficient index accepted");
    rejects<std::out_of_range>([&] { (void)view[3]; }, "coefficient index past the extent accepted");
    d::ComplexPolynomial moved(std::move(original));
    require(view[0] == moved.coefficients()[0], "moving a polynomial changed its retained coefficient storage");
    require(original.coefficients().size() == 0, "moved polynomial supplied its transferred coefficients");
    d::PolynomialSnapshot snapshot(std::move(moved));
    require(snapshot.coefficients()[0] == view[0], "owned snapshot construction copied coefficient storage");
    auto retained = std::move(snapshot);
    // NOLINTNEXTLINE(bugprone-use-after-move): the snapshot checks access after its storage was transferred.
    rejects<std::logic_error>([&] { (void)snapshot.coefficients(); }, "moved snapshot supplied a coefficient view");
    const auto borrowed = retained.coefficients();
    const auto copied_owner = retained;
    retained = d::PolynomialSnapshot(d::ComplexPolynomial{});
    require(retained.coefficients().size() == 0, "replacement snapshot kept its old coefficient extent");
    require(borrowed[0] == copied_owner.coefficients()[0] && acb_equal_si(borrowed[2], 3),
            "a retained snapshot did not keep borrowed coefficient storage alive");
}
static void buffers() {
    rejects<std::length_error>([] { d::RootVector invalid(-1); }, "negative root buffer accepted");
    rejects<std::length_error>([] { d::RealVector invalid(-1); }, "negative real buffer accepted");
    rejects<std::invalid_argument>([] { d::ConstRootView invalid(nullptr, 1); }, "null populated view accepted");
    d::RootVector roots(2);
    acb_set_si(roots[0], 7);
    acb_set_si(roots[1], 11);
    d::RootVector moved(std::move(roots));
    // NOLINTNEXTLINE(bugprone-use-after-move): the buffer explicitly promises an empty source after a move.
    require(roots.size() == 0 && roots.data() == nullptr && moved.size() == 2, "buffer transfer lost its ownership state");
    const auto &read_only = moved;
    require(acb_equal_si(read_only[0], 7) && acb_equal_si(read_only.view()[1], 11), "buffer transfer changed values");
    rejects<std::out_of_range>([&] { (void)roots[0]; }, "moved buffer accepted indexing");
    rejects<std::out_of_range>([&] { (void)moved[-1]; }, "negative index accepted");
    rejects<std::out_of_range>([&] { (void)read_only[2]; }, "const root buffer index unchecked");
    d::RealVector reals(1);
    arb_one(reals[0]);
    const auto &const_reals = reals;
    require(arb_is_one(const_reals[0]), "const real view changed value");
    rejects<std::out_of_range>([&] { (void)const_reals[1]; }, "const real buffer index unchecked");
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        buffers();
        real_polynomials();
        snapshots();
        coefficient_views();
        std::cout << "Ownership and retained-source checks passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
