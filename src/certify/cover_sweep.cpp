#include "certify/cover_sweep.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <set>

namespace certpolysolve::detail {
namespace {
int compare(size_t a, size_t b, CoverSweepWork *work) {
    if (work) ++work->comparisons;
    return (a > b) - (a < b);
}

std::vector<size_t> endpoint_ranks(const std::vector<Float> &values, CoverSweepWork *work) {
    const size_t count = values.size();
    std::vector<size_t> order(count), merged(count), ranks(count);
    std::iota(order.begin(), order.end(), size_t{0});
    const auto bits = [&](size_t i) { return std::max<slong>(1, arf_bits(values[i].value())); };
    if (work)
        for (size_t i = 0; i < count; ++i)
            work->endpoint_mantissa_bits += static_cast<uint64_t>(bits(i));
    const auto compare_values = [&](size_t a, size_t b) {
        if (work) {
            ++work->comparisons;
            ++work->rank_comparisons;
            work->rank_mantissa_bits += static_cast<uint64_t>(std::min(bits(a), bits(b)));
        }
        return arf_cmp(values[a].value(), values[b].value());
    };
    // Compare variable-length dyadics only here. A comparison examines their
    // common leading mantissa prefix; its shorter-length cost is charged to
    // the emitted item, once per merge level. Subsequent sweep queries use
    // ranks, so even repeated searches cannot repeat long arithmetic.
    for (size_t width = 1; width < count;) {
        for (size_t begin = 0; begin < count;) {
            const size_t middle = begin + std::min(width, count - begin);
            const size_t end = middle + std::min(width, count - middle);
            size_t left = begin, right = middle;
            for (size_t out = begin; out < end; ++out) {
                const bool take_left = right == end ||
                    (left < middle && compare_values(order[left], order[right]) <= 0);
                if (take_left) merged[out] = order[left++];
                else merged[out] = order[right++];
            }
            begin = end;
        }
        order.swap(merged);
        if (width >= count - width) break; // final level, without overflowing 2*width
        width *= 2;
    }
    size_t rank = 0;
    for (size_t i = 0; i < count; ++i) {
        if (i && compare_values(order[i - 1], order[i]) != 0) ++rank;
        ranks[order[i]] = rank;
    }
    return ranks;
}

// Check the exact sum's bit span with compact integer exponents BEFORE adding.
// A zero operand needs no exponent alignment. The final +1 pays a carry.
bool sum_fits(const Float &a, const Float &b, const Integer &budget) {
    if (arf_is_zero(a.value()) || arf_is_zero(b.value())) return true;
    Integer a_top, b_top, a_bottom, b_bottom, span;
    arf_abs_bound_lt_2exp_fmpz(a_top.value(), a.value());
    arf_abs_bound_lt_2exp_fmpz(b_top.value(), b.value());
    fmpz_sub_ui(a_bottom.value(), a_top.value(), static_cast<ulong>(arf_bits(a.value())));
    fmpz_sub_ui(b_bottom.value(), b_top.value(), static_cast<ulong>(arf_bits(b.value())));
    const auto top = fmpz_cmp(a_top.value(), b_top.value()) > 0 ? a_top.value() : b_top.value();
    const auto bottom = fmpz_cmp(a_bottom.value(), b_bottom.value()) < 0 ? a_bottom.value() : b_bottom.value();
    fmpz_sub(span.value(), top, bottom);
    fmpz_add_ui(span.value(), span.value(), 1);
    return fmpz_fits_si(span.value()) && fmpz_cmp(span.value(), budget.value()) <= 0;
}
} // namespace

std::optional<CoverSweep> CoverSweep::prepare(const std::vector<Root> &roots, slong precision,
                                             CoverSweepWork *work) {
    if (precision <= 0) return std::nullopt;
    // Root::publish uses p=P+16, saturated exactly as the publisher does.
    // Two distinct p-bit ends have nonzero midpoint/radius top-exponent gap
    // at most p+2. The midpoint has <=p bits, the radius <=MAG_BITS bits;
    // 2p+MAG_BITS+3 therefore includes their alignment and one addition carry.
    // Zeros are handled separately. This grows with P, not a fixed cutoff.
    Integer budget;
    fmpz_set_si(budget.value(), resources::precision_sum(precision, 16));
    fmpz_mul_ui(budget.value(), budget.value(), 2);
    fmpz_add_ui(budget.value(), budget.value(), MAG_BITS + 3);
    CoverSweep sweep(roots.size());
    // The exact values and all temporary indices belong to this call. Only
    // their order ranks survive; no Float reference escapes or crosses a move.
    std::vector<Float> real_ends(2 * roots.size()), imaginary_ends(4 * roots.size());
    for (size_t i = 0; i < roots.size(); ++i) {
        if (!acb_is_finite(roots[i].box.value())) return std::nullopt;
        for (size_t coordinate = 0; coordinate < 2; ++coordinate) {
            const auto part = coordinate ? acb_imagref(roots[i].box.value()) : acb_realref(roots[i].box.value());
            Float midpoint, radius;
            arf_set(midpoint.value(), arb_midref(part));
            arf_set_mag(radius.value(), arb_radref(part));
            if (!sum_fits(midpoint, radius, budget)) return std::nullopt;
            auto &ends = coordinate ? imaginary_ends : real_ends;
            const size_t offset = coordinate ? 4 * i : 2 * i;
            arf_sub(ends[offset].value(), midpoint.value(), radius.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_add(ends[offset + 1].value(), midpoint.value(), radius.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
        }
        arf_neg(imaginary_ends[4 * i + 2].value(), imaginary_ends[4 * i + 1].value());
        arf_neg(imaginary_ends[4 * i + 3].value(), imaginary_ends[4 * i].value());
    }
    const auto real_ranks = endpoint_ranks(real_ends, work);
    const auto imaginary_ranks = endpoint_ranks(imaginary_ends, work);
    for (size_t i = 0; i < roots.size(); ++i) {
        sweep.rectangles_[i].x = {real_ranks[2 * i], real_ranks[2 * i + 1]};
        for (size_t family = 0; family < 2; ++family)
            sweep.rectangles_[i].y[family] = {imaginary_ranks[4 * i + 2 * family],
                                             imaginary_ranks[4 * i + 2 * family + 1]};
    }
    return sweep;
}

struct CoverSweep::Event {
    size_t id;
    size_t family;
    bool start;
};

std::vector<CoverSweep::Event> CoverSweep::events(bool reflected, CoverSweepWork *work) const {
    std::vector<Event> out;
    for (size_t i = 0; i < rectangles_.size(); ++i)
        for (size_t family = 0; family < (reflected ? 2u : 1u); ++family) {
            out.push_back({i, family, true});
            out.push_back({i, family, false});
        }
    std::sort(out.begin(), out.end(), [&](const Event &a, const Event &b) {
        const int cmp = compare(rectangles_[a.id].x[a.start ? 0 : 1],
                                rectangles_[b.id].x[b.start ? 0 : 1], work);
        if (cmp) return cmp < 0;
        if (a.start != b.start) return a.start; // closed: every start precedes every end
        if (a.id != b.id) return a.id < b.id;
        return a.family < b.family;
    });
    return out;
}

struct CoverSweep::IntervalLess {
    // These borrows stay inside a sweep call; the owned rectangles never move
    // while a set or a query refers to them. No view escapes to solver state.
    const CoverSweep &geometry;
    size_t family;
    CoverSweepWork *work;
    struct Start { size_t lower; };
    using is_transparent = void;
    const std::array<size_t, 2> &interval(size_t id) const { return geometry.rectangles_[id].y[family]; }
    bool operator()(size_t a, size_t b) const {
        const int cmp = compare(interval(a)[0], interval(b)[0], work);
        return cmp ? cmp < 0 : a < b;
    }
    bool operator()(size_t id, const Start &x) const { return compare(interval(id)[0], x.lower, work) < 0; }
    bool operator()(const Start &x, size_t id) const { return compare(x.lower, interval(id)[0], work) < 0; }

    // Active intervals in this family are disjoint. Only the predecessor
    // and first following interval can meet, plus one extra if self is skipped.
    template <class Set>
    auto overlap(const Set &active, const std::array<size_t, 2> &query,
                 size_t exclude = std::numeric_limits<size_t>::max()) const {
        auto next = active.lower_bound(Start{query[0]});
        if (next != active.begin()) {
            const auto previous = std::prev(next);
            if (*previous != exclude && compare(interval(*previous)[1], query[0], work) >= 0)
                return previous;
        }
        for (; next != active.end() && compare(interval(*next)[0], query[1], work) <= 0; ++next)
            if (*next != exclude) return next;
        return active.end();
    }
};

bool CoverSweep::disjoint(CoverSweepWork *work) const {
    const IntervalLess less{*this, 0, work};
    std::set<size_t, IntervalLess> active(less);
    for (const auto &event : events(false, work)) {
        if (!event.start) active.erase(event.id);
        else {
            if (less.overlap(active, rectangles_[event.id].y[0]) != active.end()) return false;
            active.insert(event.id);
        }
    }
    return true;
}

std::vector<bool> CoverSweep::reflection_met(CoverSweepWork *work) const {
    const std::array<IntervalLess, 2> less{{{*this, 0, work}, {*this, 1, work}}};
    using Set = std::set<size_t, IntervalLess>;
    std::array<Set, 2> all{Set(less[0]), Set(less[1])};
    std::array<Set, 2> unmarked{Set(less[0]), Set(less[1])};
    std::vector<bool> met(rectangles_.size(), false);
    const auto mark = [&](size_t id) {
        if (met[id]) return;
        met[id] = true;
        if (work) ++work->marked;
        for (auto &family : unmarked) family.erase(id);
    };
    for (const auto &event : events(true, work)) {
        const size_t id = event.id, family = event.family, opposite = 1 - family;
        if (!event.start) {
            all[family].erase(id);
            unmarked[family].erase(id);
            continue;
        }
        const auto &query = rectangles_[id].y[family];
        if (less[opposite].overlap(all[opposite], query, id) != all[opposite].end()) mark(id);
        while (true) {
            const auto found = less[opposite].overlap(unmarked[opposite], query, id);
            if (found == unmarked[opposite].end()) break;
            mark(*found); // removes this witness once, even if there are quadratically many cross edges
        }
        all[family].insert(id);
        if (!met[id]) unmarked[family].insert(id);
    }
    return met;
}
} // namespace certpolysolve::detail
