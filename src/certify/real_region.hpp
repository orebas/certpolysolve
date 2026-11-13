#pragma once
// D1, Part 2 (docs/real-interval-contract-note-20261015.md): the real roots of a real-only query in the closed
// region of cps_options_set_region. The region's trace on the real axis is the closed interval [lower, upper]
// (either side may be unbounded) when the region meets the axis, and empty otherwise. Every answer of the real-only
// routes (real_only_route.cpp) holds each real root alone in its exported interval; a root is kept when that
// interval lies in the closed interval, dropped when it lies outside, and decided exactly when the interval holds an
// end: the end is the root when the input vanishes there (the exact test below), and otherwise the root lies
// strictly on one side. For a root of odd multiplicity the input changes sign once in the interval, at the root, so
// its signs at the end and at the interval's lower end tell the side, and the interval is cut at a dyadic point
// between the end and the root (place, cut). A root of even multiplicity has no sign change to read: it stays
// undecided, and the route solves again at doubled accuracy until its interval lies on one side. History:
// docs/code-history-20261021.md, real_region.hpp 1.
#include "core/cover.hpp"
#include "core/query.hpp"
#include "certify/real_sign_refinement.hpp"
#include "numeric/evaluate/sparse_real_evaluator.hpp"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <optional>
#include <utility>
#include <vector>

namespace certpolysolve::detail {

class RealRegionFilter {
  public:
    // `o`'s region and exact storage policy; `terms` are the input's nonzero terms by ascending exponent (exact,
    // real). The ends become exact rationals here, under the policy of the input's own exact import: refusal()
    // names an end it refuses.
    RealRegionFilter(const cps_options &o, std::vector<SparseRealTerm> terms) : terms_(std::move(terms)) {
        meets_axis_ = (!o.region.imag_lower || o.region.imag_lower->sign() <= 0) &&
                      (!o.region.imag_upper || o.region.imag_upper->sign() >= 0);
        refusal_ = exact_ends(o, ends_);
    }
    // The region's real ends as exact rationals, under the options' exact storage policy (one allowance for both
    // ends, the lower first), into `ends`; the end it refuses, named, or nullptr. The options alone decide it, so
    // the solve checks it before any work (C8, c_api.cpp) and a refused end is a usage error there.
    static const char *exact_ends(const cps_options &o, std::array<std::optional<Rational>, 2> &ends) {
        uint64_t budget = o.exact_storage_limit();
        const char *refused = nullptr;
        const auto end = [&](const std::optional<Number> &given, std::optional<Rational> &out, const char *name) {
            if (!given) return;
            Rational value;
            if (given->rational_with_budget(value.value(), budget, o.allow_exponent_expansion)) out = std::move(value);
            else refused = name;
        };
        end(o.region.real_lower, ends[0], "the region's real lower end exceeds the exact storage policy "
                                          "(the resource limits' exact_materialization_bits and "
                                          "CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION allow its expansion)");
        end(o.region.real_upper, ends[1], "the region's real upper end exceeds the exact storage policy "
                                          "(the resource limits' exact_materialization_bits and "
                                          "CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION allow its expansion)");
        return refused;
    }
    const char *refusal() const { return refusal_; }
    bool meets_axis() const { return meets_axis_; }

    struct Outcome {
        size_t undecided = 0;         // roots kept whose interval holds an end that is not their root
        const char *refusal = nullptr; // an exact test the exact storage policy refuses
    };
    // Keeps the roots of a real-only answer that lie in the interval, and the undecided ones, so that the answer at
    // a stop holds every root that may lie in it. An individual root (its multiplicity established: the only kind
    // a goal-met answer holds) is decided exactly against an end it reaches; one equal to a dyadic end becomes
    // that exact point, one equal to another end keeps its interval, which holds it alone; one of odd multiplicity
    // whose interval holds an end it is not equal to is placed by the input's signs and its interval cut to the
    // region's side, or to the root's point where the input vanishes at the cut. A counted cluster (an answer at a
    // stop) is only dropped when its interval lies wholly outside, and is otherwise kept as it is, undecided: the
    // exact test says nothing about its other roots. `precision` is the answer's working precision, where the
    // signs' evaluations start. History: docs/code-history-20261021.md, real_region.hpp 2.
    Outcome keep(std::vector<Root> &roots, const cps_options &o, slong precision) {
        Outcome outcome;
        std::vector<Root> kept;
        for (auto &root : roots) {
            if (!meets_axis_) continue;
            const bool individual = root.multiplicity == root.count;
            // Compare the compact dyadics directly. A distant root or a stopped cluster needs no expanded
            // rational endpoints; expanding 1@E here would allocate in proportion to |E| before any exact test.
            Verdict position = Verdict::Inside;
            bool touches = false;
            for (size_t which = 0; which < ends_.size(); ++which) {
                if (!ends_[which]) continue;
                const int lo_e = compare(root.endpoints[0], *ends_[which]);
                const int hi_e = compare(root.endpoints[1], *ends_[which]);
                const bool lower = which == 0;
                if (lower ? hi_e < 0 : lo_e > 0) {
                    position = Verdict::Outside;
                    break;
                }
                if (lower ? lo_e < 0 : hi_e > 0) position = Verdict::Undecided;
                touches = touches || (lo_e <= 0 && hi_e >= 0);
            }
            if (position == Verdict::Outside) continue;
            if (!individual || (position == Verdict::Inside && !touches)) {
                if (position == Verdict::Undecided) ++outcome.undecided;
                kept.push_back(std::move(root));
                continue;
            }
            uint64_t remaining = o.exact_storage_limit();
            auto lo_exact = exported(root.endpoints[0], o, remaining);
            auto hi_exact = exported(root.endpoints[1], o, remaining);
            if (!lo_exact || !hi_exact) {
                // A wholly contained interval already proves membership, even if an optional exact-point
                // reduction is too large. An interval crossing an end must stay undecided at this resource stop.
                if (position == Verdict::Undecided) {
                    ++outcome.undecided;
                    outcome.refusal = "placing a root against the region's end needs expanded exact endpoints "
                                      "which the exact storage or compact-exponent policy refuses";
                }
                kept.push_back(std::move(root));
                continue;
            }
            const Rational lo_given = std::move(*lo_exact), hi_given = std::move(*hi_exact);
            Rational lo = lo_given, hi = hi_given; // cut by place against an end the interval holds
            Verdict verdict = Verdict::Inside;
            const Rational *at = nullptr;
            for (size_t which = 0; which < ends_.size() && verdict != Verdict::Outside; ++which) {
                if (!ends_[which]) continue;
                bool equal = false;
                auto side = against(lo, hi, which, root.multiplicity, equal, o, precision);
                if (side == Verdict::Refused) {
                    // The root stays, undecided: the answer at this stop holds every root that may lie in the region.
                    outcome.refusal = "an end of the region lies in a root's interval, and its exact test needs the "
                                      "input in exact integer form, which the exact storage policy refuses "
                                      "(the resource limits' exact_materialization_bits)";
                    side = Verdict::Undecided;
                }
                if (equal) at = &*ends_[which];
                if (side == Verdict::Outside || verdict == Verdict::Inside) verdict = side;
            }
            if (verdict == Verdict::Outside) continue;
            if (verdict == Verdict::Undecided) ++outcome.undecided;
            const bool moved = !fmpq_equal(lo.value(), lo_given.value()) || !fmpq_equal(hi.value(), hi_given.value());
            if (at && dyadic(*at)) publish_point(root, *at, "; the region's end, where the input vanishes");
            else if (moved && fmpq_equal(lo.value(), hi.value()))
                publish_point(root, lo, "; placed against the region's end by the input's signs, which vanish here");
            else if (moved)
                publish_cut(root, lo, hi);
            kept.push_back(std::move(root));
        }
        roots = std::move(kept);
        return outcome;
    }

  private:
    enum class Verdict { Inside, Outside, Undecided, Refused };
    enum class Vanishes { Yes, No, Refused, Unknown };
    std::vector<SparseRealTerm> terms_;
    std::array<std::optional<Rational>, 2> ends_; // lower, upper
    bool meets_axis_ = true;
    const char *refusal_ = nullptr;
    // The exact test's state: the input as an integer polynomial along its support (built at the first test), and
    // each end's verdict once known.
    bool imported_ = false;
    std::vector<std::pair<slong, Integer>> integer_;
    std::array<std::optional<bool>, 2> vanishes_;
    // The signs' evaluator over the input's terms (built at the first sign).
    std::optional<SparseRealEvaluator> evaluator_;

    // Compare x with r/s by comparing the exact binary float s*x with the integer r. Its significand grows
    // with s, already admitted as a region end; x's exponent stays compact regardless of its magnitude.
    static int compare(const Endpoint &x, const Rational &end) {
        Float scaled, numerator;
        arf_mul_fmpz(scaled.value(), x.value(), fmpq_denref(end.value()), ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_set_fmpz(numerator.value(), fmpq_numref(end.value()));
        return arf_cmp(scaled.value(), numerator.value());
    }
    // Only a boundary-placement test needs rational endpoints. Apply the same exact-storage and compact-expansion
    // policy as coefficient preparation, with one shared storage allowance for both endpoints.
    static std::optional<Rational> exported(const Endpoint &x, const cps_options &o, uint64_t &remaining) {
        Number compact;
        Integer mantissa;
        x.get(mantissa.value(), compact.exponent.value());
        fmpq_set_fmpz(compact.significand.value(), mantissa.value());
        const bool expand = o.allow_exponent_expansion ||
                            expansion_fits(compact.owned_storage_bits(), compact.expanded_storage_bits(),
                                           compact.zero() ? 0 : 1, o.exact_expansion_allowance());
        Rational out;
        if (!compact.rational_with_budget(out.value(), remaining, expand)) return std::nullopt;
        return out;
    }
    // A power of two in the denominator (it is positive): its lowest set bit is its only one.
    static bool dyadic(const Rational &x) {
        const fmpz *denominator = fmpq_denref(x.value());
        return fmpz_val2(denominator) + 1 == fmpz_bits(denominator);
    }
    // The precision at which a dyadic x is held and exported exactly (its bits and a guard).
    static slong exact_precision(const Rational &x) {
        const slong bits = static_cast<slong>(fmpz_bits(fmpq_numref(x.value())) + fmpz_bits(fmpq_denref(x.value())));
        return std::max<slong>(64, resources::precision_sum(bits, 16));
    }
    // The root is the dyadic point x (an end of the region, or a cut point): its box becomes the exact point.
    static void publish_point(Root &root, const Rational &x, const char *note) {
        const slong precision = exact_precision(x);
        arb_set_fmpq(acb_realref(root.box.value()), x.value(), precision);
        arb_zero(acb_imagref(root.box.value()));
        if (!arb_is_exact(acb_realref(root.box.value()))) throw std::logic_error("a dyadic end is not exact");
        root.realness = CPS_REAL;
        root.certificate += note;
        root.publish(precision);
    }
    // A dyadic x = numerator 2^-k as an exact binary float.
    static void dyadic_float(Float &out, const Rational &x) {
        Integer exponent;
        fmpz_set_si(exponent.value(), 1 - static_cast<slong>(fmpz_bits(fmpq_denref(x.value()))));
        arf_set_fmpz_2exp(out.value(), fmpq_numref(x.value()), exponent.value());
    }
    // The root's interval cut to [lo, hi] (dyadic ends): exported as these very ends, its box a ball that holds them,
    // as Root::publish leaves it (the box holds the exported box; publish itself would round the ends outward, past
    // the region's end that one of them may be).
    static void publish_cut(Root &root, const Rational &lo, const Rational &hi) {
        Float a, b;
        dyadic_float(a, lo);
        dyadic_float(b, hi);
        root.endpoints = {Endpoint(a.value()), Endpoint(b.value()), Endpoint{}, Endpoint{}};
        arb_set_interval_arf(acb_realref(root.box.value()), a.value(), b.value(),
                             std::max(exact_precision(lo), exact_precision(hi)));
        arb_zero(acb_imagref(root.box.value()));
        root.realness = CPS_REAL;
        root.certificate += "; placed against the region's end by the input's signs";
    }
    // Where the root alone in [lo, hi], of `multiplicity`, lies against end `which` (0 the lower end: the root must be
    // at least it; 1 the upper: at most it). `equal` is set when the root is the end; place may cut [lo, hi].
    Verdict against(Rational &lo, Rational &hi, size_t which, uint64_t multiplicity, bool &equal, const cps_options &o,
                    slong precision) {
        const Rational &e = *ends_[which];
        const bool lower = which == 0;
        const int lo_e = fmpq_cmp(lo.value(), e.value()), hi_e = fmpq_cmp(hi.value(), e.value());
        if (lower ? hi_e < 0 : lo_e > 0) return Verdict::Outside; // wholly beyond the end
        if (lower ? lo_e > 0 : hi_e < 0) return Verdict::Inside;  // wholly on the region's side
        if (lo_e == 0 && hi_e == 0) {                              // the interval is the end's point
            equal = true;
            return Verdict::Inside;
        }
        // The interval holds the end. Where the end is the interval's own end on the region's side ([e, hi] for a
        // lower end, [lo, e] for an upper one) the root is in the region whatever the exact test says; the test
        // only shows whether the root is the end, to be published as its point.
        const bool touching = lower ? lo_e == 0 : hi_e == 0;
        switch (vanishes(which, o)) {
        case Vanishes::Yes: equal = true; return Verdict::Inside; // the interval's one root is the end
        case Vanishes::Refused: return touching ? Verdict::Inside : Verdict::Refused;
        case Vanishes::Unknown: return touching ? Verdict::Inside : Verdict::Undecided; // cancelled: the session stops
        case Vanishes::No: break;
        }
        if (touching) return Verdict::Inside;
        // Not the end: an interval that reaches the end only from outside, [lo, e] below a lower end or [e, hi] above
        // an upper one, holds its root outside; otherwise the root lies strictly on one side of the end, which the
        // input's signs tell where it changes sign at the root (odd multiplicity), and refinement otherwise.
        if (lower ? hi_e == 0 : lo_e == 0) return Verdict::Outside;
        if (multiplicity % 2 == 1) return place(lo, hi, which, o, precision);
        return Verdict::Undecided;
    }
    // The root alone in [lo, hi], of odd multiplicity, against end e strictly inside, where the input does not vanish:
    // the input changes sign once in [lo, hi], at the root, so the root lies in (lo, e) exactly when the signs at lo
    // and at e differ. A root on the region's side gets its interval cut there (cut); one at lo itself (the input
    // vanishes there) becomes that point. Undecided where a sign is (cancelled, or the exact value refused).
    Verdict place(Rational &lo, Rational &hi, size_t which, const cps_options &o, slong precision) {
        const Rational &e = *ends_[which];
        const bool lower = which == 0;
        bool zero = false;
        const int at_e = sign_at(e, o, precision, zero);
        if (zero) throw std::logic_error("the input vanishes at an end its exact test found not to be a root");
        if (!at_e) return Verdict::Undecided;
        const int at_lo = sign_at(lo, o, precision, zero);
        if (zero) { // the root is lo, below e
            hi = lo;
            return lower ? Verdict::Outside : Verdict::Inside;
        }
        if (!at_lo) return Verdict::Undecided;
        const bool below = at_lo != at_e;
        if (lower == below) return Verdict::Outside; // below a lower end, or above an upper one
        return cut(lo, hi, which, at_e, o, precision) ? Verdict::Inside : Verdict::Undecided;
    }
    // The root lies on the region's side of e (above a lower end, below an upper one), strictly. Moves the interval's
    // end on the other side of e (lo for a lower end, hi for an upper one) to a dyadic point c between e and the root:
    // c is e rounded toward the root to a multiple of 2^-k, so it lies within 2^-k of e, and k starts where that
    // multiple lies inside the interval (k0: 2^-k0 below the distance from e to the far end, the interval's other
    // end). Where the input's sign at c is its sign at e (`at_e`), no root lies between them, and the interval runs
    // from c to the far end; where the input vanishes at c, the root is c; otherwise the root lies between e and c, c
    // becomes the far end, and k goes twice as many bits below k0, plus one (0, 1, 3, 7, ... bits: a root 2^-d from
    // e takes about log2(d) steps). It ends because the root is not e. The search stops undecided where a cut point
    // would need more bits than the exact storage policy allows (the region's ends obey it too), and where a sign is
    // undecided (cancelled, or the exact value refused). The precision cap bounds the signs' ball evaluations
    // (sign_at), not the search: beyond it the exact value decides, as the exact test at an end does.
    bool cut(Rational &lo, Rational &hi, size_t which, int at_e, const cps_options &o, slong precision) {
        const Rational &e = *ends_[which];
        const bool lower = which == 0;
        Rational far = lower ? hi : lo, distance, scaled, c;
        fmpq_sub(distance.value(), far.value(), e.value());
        fmpq_abs(distance.value(), distance.value());
        // distance > 2^(t - 1) with t the difference of its numerator's and denominator's bits, so 2^-k0 < distance
        // at k0 = 2 - t
        const slong k0 = 2 - (static_cast<slong>(fmpz_bits(fmpq_numref(distance.value()))) -
                              static_cast<slong>(fmpz_bits(fmpq_denref(distance.value()))));
        const uint64_t reach = o.exact_storage_limit();
        Integer rounded;
        // (A grid finer than a precision can name is beyond any search: the root's distance from e is a positive
        // number of the input's and e's own sizes, so the search ends far sooner.)
        const slong finest = static_cast<slong>(resources::precision_max / 2);
        for (slong below = 0;; below = resources::precision_sum(resources::precision_product(below, 2), 1)) {
            if (static_cast<uint64_t>(below) > reach || below >= finest - std::abs(k0) || o.cancelled()) return false;
            const slong k = k0 + below;
            // c = e rounded toward the root (up for a lower end, down for an upper one) to a multiple of 2^-k
            if (k >= 0) fmpq_mul_2exp(scaled.value(), e.value(), static_cast<ulong>(k));
            else fmpq_div_2exp(scaled.value(), e.value(), static_cast<ulong>(-k));
            if (lower) fmpz_cdiv_q(rounded.value(), fmpq_numref(scaled.value()), fmpq_denref(scaled.value()));
            else fmpz_fdiv_q(rounded.value(), fmpq_numref(scaled.value()), fmpq_denref(scaled.value()));
            fmpq_set_fmpz(c.value(), rounded.value());
            if (k >= 0) fmpq_div_2exp(c.value(), c.value(), static_cast<ulong>(k));
            else fmpq_mul_2exp(c.value(), c.value(), static_cast<ulong>(-k));
            if (fmpq_equal(c.value(), far.value())) continue; // the finer grid's point is the far end's: a finer one
            bool zero = false;
            // c = e (e is a multiple of 2^-k): the interval reaches e, which is not its root.
            const int at_c = fmpq_equal(c.value(), e.value()) ? at_e : sign_at(c, o, precision, zero);
            if (zero) {
                lo = c;
                hi = c;
                return true;
            }
            if (!at_c) return false;
            if (at_c == at_e) {
                if (lower) lo = c, hi = far;
                else lo = far, hi = c;
                return true;
            }
            far = c;
        }
    }
    // The input's sign at the exact point x, from its exact form (imported under the exact storage policy; undecided
    // without it). A point the rational root theorem allows to be a root gets its exact value (exact_sign), since no
    // ball shows a zero; any other point gets certified evaluations along the input's terms (SparseRealEvaluator, as
    // the axis charts' and the sparse recursion's signs), x held at the evaluation's precision (exactly when dyadic),
    // from the answer's working precision doubling up to the query's cap, then the exact value. `zero` is set where
    // the input vanishes at x; 0 without it is undecided (cancelled, or the exact form refused).
    int sign_at(const Rational &x, const cps_options &o, slong precision, bool &zero) {
        zero = false;
        if (!import(o)) return 0;
        if (may_vanish(x)) {
            const auto exact = exact_sign(x, o);
            if (!exact) return 0;
            zero = *exact == 0;
            return *exact;
        }
        if (!evaluator_) evaluator_.emplace(terms_);
        const slong cap = std::max(real_sign_cap(o), precision);
        RealBall point, value, slope;
        for (slong bits = precision;; bits = std::min(resources::precision_product(bits, 2), cap)) {
            if (o.cancelled()) return 0;
            arb_set_fmpq(point.value(), x.value(), std::max(bits, exact_precision(x)));
            evaluator_->evaluate2(value.value(), slope.value(), point.value(), bits);
            if (const int sign = real_sign_of(value.value())) return sign;
            if (bits >= cap) break;
        }
        const auto exact = exact_sign(x, o);
        return exact ? *exact : 0;
    }
    // The rational root theorem on the imported integer form P: x = r/s in lowest terms is a root only if s divides
    // P's leading coefficient and r its lowest nonzero one; 0 is a root exactly when the lowest exponent is positive.
    bool may_vanish(const Rational &x) const {
        if (fmpq_is_zero(x.value())) return integer_.front().first > 0;
        return fmpz_divisible(integer_.back().second.value(), fmpq_denref(x.value())) &&
               fmpz_divisible(integer_.front().second.value(), fmpq_numref(x.value()));
    }
    // The input's exact sign at x, from its integer form (scaled_value); nothing where the exact storage policy
    // refuses that form or the query is cancelled.
    std::optional<int> exact_sign(const Rational &x, const cps_options &o) {
        if (terms_.empty() || !import(o)) return std::nullopt;
        // At 0 the lowest term decides: the constant's sign, or zero (the zero root).
        if (fmpq_is_zero(x.value()))
            return integer_.front().first > 0 ? 0 : fmpz_sgn(integer_.front().second.value());
        Integer h;
        if (!scaled_value(x, h, o)) return std::nullopt;
        // h = s^(k_m - k_0) x^(-k_0) P(x) with s > 0, and P = D p with D > 0: p(x) has h's sign times x's to the
        // power k_0.
        const int sign = fmpz_sgn(h.value());
        return integer_.front().first % 2 == 1 && fmpq_sgn(x.value()) < 0 ? -sign : sign;
    }
    // Whether the input vanishes at end `which`, exactly.
    Vanishes vanishes(size_t which, const cps_options &o) {
        if (vanishes_[which]) return *vanishes_[which] ? Vanishes::Yes : Vanishes::No;
        const Rational &e = *ends_[which];
        if (terms_.empty()) return Vanishes::Unknown; // the zero polynomial: the routes stop before any root
        if (fmpq_is_zero(e.value())) {
            vanishes_[which] = terms_.front().exponent > 0; // the zero root
            return *vanishes_[which] ? Vanishes::Yes : Vanishes::No;
        }
        if (!import(o)) return Vanishes::Refused;
        if (!may_vanish(e)) { // the rational root theorem
            vanishes_[which] = false;
            return Vanishes::No;
        }
        Integer h;
        if (!scaled_value(e, h, o)) return Vanishes::Unknown;
        vanishes_[which] = fmpz_is_zero(h.value());
        return *vanishes_[which] ? Vanishes::Yes : Vanishes::No;
    }
    // For x = r/s (nonzero, lowest terms) and the imported integer form P: Horner's rule along the support from the
    // top, with integers: H = P_(k_m), then for each lower exponent H <- r^(k_(i+1) - k_i) H + s^(k_m - k_i) P_(k_i),
    // which ends at s^(k_m - k_0) P(r/s) / (r/s)^(k_0): zero exactly when P(r/s) is (r is not zero). One step per
    // nonzero term, the powers along the gaps. False when cancelled.
    bool scaled_value(const Rational &x, Integer &h, const cps_options &o) const {
        const fmpz *r = fmpq_numref(x.value()), *s = fmpq_denref(x.value());
        Integer scale, power;
        fmpz_set(h.value(), integer_.back().second.value());
        fmpz_one(scale.value());
        for (size_t i = integer_.size() - 1; i-- > 0;) {
            if (o.cancelled()) return false;
            const ulong gap = static_cast<ulong>(integer_[i + 1].first - integer_[i].first);
            fmpz_pow_ui(power.value(), r, gap);
            fmpz_mul(h.value(), h.value(), power.value());
            fmpz_pow_ui(power.value(), s, gap);
            fmpz_mul(scale.value(), scale.value(), power.value());
            fmpz_addmul(h.value(), scale.value(), integer_[i].second.value());
        }
        return true;
    }
    // The input as an integer polynomial along its support (a common denominator cleared), under the policy of the
    // solver's own exact import (streamed_real_polynomial): the options' exact storage limit,
    // and expansion of compact exponents only where they allow it or the expanded form exceeds the compact one by at
    // most the compact-exponent rule's allowance a number (expansion_fits). False where the policy refuses.
    bool import(const cps_options &o) {
        if (imported_) return !integer_.empty();
        imported_ = true;
        uint64_t owned = 0, expanded = 0, numbers = 0;
        for (const auto &term : terms_) {
            owned = resources::add(owned, term.coefficient.owned_storage_bits());
            expanded = resources::add(expanded, term.coefficient.expanded_storage_bits());
            numbers += !term.coefficient.zero();
        }
        const bool expansion =
            o.allow_exponent_expansion || expansion_fits(owned, expanded, numbers, o.exact_expansion_allowance());
        uint64_t remaining = o.exact_storage_limit();
        std::vector<Rational> rationals(terms_.size());
        Integer denominator;
        fmpz_one(denominator.value());
        for (size_t i = 0; i < terms_.size(); ++i) {
            if (!terms_[i].coefficient.rational_with_budget(rationals[i].value(), remaining, expansion)) return false;
            fmpz_lcm(denominator.value(), denominator.value(), fmpq_denref(rationals[i].value()));
        }
        // Clearing the denominators can multiply the storage by the number of terms: checked before it is done.
        if (static_cast<uint64_t>(fmpz_bits(denominator.value())) > remaining / terms_.size()) return false;
        std::vector<std::pair<slong, Integer>> integer(terms_.size());
        Integer multiplier;
        for (size_t i = 0; i < terms_.size(); ++i) {
            fmpz_divexact(multiplier.value(), denominator.value(), fmpq_denref(rationals[i].value()));
            integer[i].first = terms_[i].exponent;
            fmpz_mul(integer[i].second.value(), multiplier.value(), fmpq_numref(rationals[i].value()));
        }
        integer_ = std::move(integer);
        return true;
    }
};

} // namespace certpolysolve::detail
