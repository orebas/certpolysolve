#pragma once
#include <array>
#include <complex>
#include <limits>
#include <span>
#include <vector>

namespace certpolysolve::detail {
// A node is a far field for a point whose distance from the node's centre is
// at least its radius over this ratio.
inline constexpr double interaction_acceptance = 0.20;
// The fewest series terms whose truncation bound q^k / (1 - q) at the
// acceptance ratio q, relative to the node's own sum of magnitudes, is at most
// `epsilon`.
constexpr size_t interaction_series_order(double epsilon) {
    double bound = 1 / (1 - interaction_acceptance);
    size_t terms = 0;
    for (; bound > epsilon; ++terms)
        bound *= interaction_acceptance;
    return terms;
}

// The truncation bound q^order / (1 - q) of an Order-term far field at the
// acceptance ratio q, relative to the node's own sum of magnitudes (the bound
// interaction_series_order sizes the series by).
constexpr double interaction_truncation(size_t order, double acceptance = interaction_acceptance) {
    double bound = 1 / (1 - acceptance);
    for (size_t k = 0; k < order; ++k)
        bound *= acceptance;
    return bound;
}

// A tree's far-field rule, as a type, so that its acceptance ratio is a
// constant in the traversal's code: a type with `static constexpr double
// acceptance`. The coefficient routes' trees take interaction_acceptance; the
// program route's binary64 stage takes a coarser one
// (heuristics::ProgramFarField).
struct DefaultFarField {
    static constexpr double acceptance = interaction_acceptance;
};

// Candidate-only Cauchy interaction tree. Far fields use a finite multipole
// series of Order terms, taken where FarField's rule accepts them; neighbours
// and self exclusion are returned to the caller. Nothing produced by this class
// is a root certificate.
template <size_t Order, class FarField = DefaultFarField> class InteractionTreeOf {
  public:
    // Terms of each far field's multipole series.
    static constexpr size_t order = Order;

  private:
    struct Node {
        size_t begin = 0, end = 0, left = 0, right = 0;
        std::complex<double> center{};
        // The centre's second part, with lows_ (the centre is center + center_low).
        std::complex<double> center_low{};
        double radius = 0;
        // The largest rounding error among its points (errors_).
        double error = 0;
        std::array<std::complex<double>, order> moments{};
    };
    std::vector<std::complex<double>> points_;
    // When the caller's points are not binary64 (the multiprecision packets,
    // prepare_interactions): each point's second double, what its first one
    // leaves of it (a point is points_[j] + lows_[j]), and a bound on what the
    // two leave, |dx| + |dy|. Empty: the copies are the points.
    std::vector<std::complex<double>> lows_;
    std::vector<double> errors_;
    // A point's displacement from a centre (high and low parts), each part's
    // difference taken first: the high parts of nearby points differ exactly.
    std::complex<double> displacement(size_t point, std::complex<double> center,
                                      std::complex<double> center_low) const {
        return (points_[point] - center) + (lows_[point] - center_low);
    }
    std::vector<size_t> indices_, positions_;
    std::vector<Node> nodes_;
    size_t build_terms_ = 0;
    size_t build(size_t begin, size_t end);
    size_t build_two_part(size_t index, size_t begin, size_t end);
    bool visit(size_t node, size_t target, std::complex<double>& far, std::vector<size_t>& near,
               double& tail, size_t& far_fields) const;

  public:
    explicit InteractionTreeOf(std::span<const std::complex<double>> points, std::span<const double> errors = {},
                               std::span<const std::complex<double>> lows = {});
    bool evaluate(size_t target, std::complex<double>& far, std::vector<size_t>& near, double& tail) const;
    // The same, counting the far fields summed (each a series of `order` terms).
    bool evaluate(size_t target, std::complex<double>& far, std::vector<size_t>& near, double& tail,
                  size_t& far_fields) const;
    // The multipole terms the construction summed: every point's, in every
    // node that holds it.
    size_t build_terms() const { return build_terms_; }
    // Moves a point. The moments of every node that holds it follow, so later
    // far fields see the point where it is now; a node whose radius the point
    // leaves takes the larger radius, its moments rescaled. Returns the
    // multipole terms updated. For binary64 points (no rounding errors, no
    // second parts).
    size_t move(size_t point, std::complex<double> to);
};
// The multiprecision candidates' tree (candidates.cpp).
using InteractionTree = InteractionTreeOf<20>;
// The binary64 Aberth packet's tree (machine_candidates.cpp): its far fields
// are as accurate as the direct sum they replace, truncated below binary64's
// unit roundoff (23 terms at the acceptance ratio 0.2; 20 left about 1e-14).
using MachineInteractionTree =
    InteractionTreeOf<interaction_series_order(std::numeric_limits<double>::epsilon() / 2)>;
} // namespace certpolysolve::detail
