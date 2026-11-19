#include "certify/attempts/conjugate_pairs.hpp"
#include <algorithm>
#include <cmath>

namespace certpolysolve::detail {
namespace {
// One side of the axis: its points (index into the caller's set, real part, size of the imaginary part) and the two
// orders the other side searches it by.
struct Side {
    std::vector<std::size_t> index;
    std::vector<double> re, im;
    std::vector<std::size_t> by_re, by_im; // positions in `index`, ascending in re and in im
    std::size_t size() const { return index.size(); }
    void add(std::size_t i, double real, double imaginary_size) {
        index.push_back(i);
        re.push_back(real);
        im.push_back(imaginary_size);
    }
    void sort() {
        by_re.resize(size());
        by_im.resize(size());
        for (std::size_t k = 0; k < size(); ++k)
            by_re[k] = by_im[k] = k;
        std::sort(by_re.begin(), by_re.end(), [&](std::size_t a, std::size_t b) { return re[a] < re[b]; });
        std::sort(by_im.begin(), by_im.end(), [&](std::size_t a, std::size_t b) { return im[a] < im[b]; });
    }
};
// A point's candidates on the other side: how many, and one of them. Not looked for (`searched` false) where the
// point has more than the scan limit of neighbours in both coordinates.
struct Candidates {
    bool searched = true;
    std::size_t count = 0, one = 0;
};
// Whether (re, im) and the other side's point k mirror each other: the test of the header, symmetric in the two.
bool mirrors(double re, double im, const Side &other, std::size_t k, int bits) {
    const double size = std::max(std::max(std::abs(re), im), std::max(std::abs(other.re[k]), other.im[k]));
    const double tolerance = std::ldexp(size, -bits);
    return std::abs(re - other.re[k]) <= tolerance && std::abs(im - other.im[k]) <= tolerance;
}
// The candidates of the point (re, im) among `other`, by a scan of all of it where it is small (`sorted` false) and
// otherwise of the shorter of its two windows: the points whose real part, or whose imaginary size, lies within
// twice the point's own tolerance of the point's. Every candidate is in both windows: its size is within a
// tolerance of the point's, so the tolerance the two share is below twice the point's own (bits >= 2).
Candidates candidates_of(double re, double im, const Side &other, bool sorted, int bits, std::size_t scan_limit) {
    Candidates out;
    const auto consider = [&](std::size_t k) {
        if (mirrors(re, im, other, k, bits)) {
            ++out.count;
            out.one = k;
        }
    };
    if (!sorted) {
        for (std::size_t k = 0; k < other.size(); ++k)
            consider(k);
        return out;
    }
    const double window = std::ldexp(std::max(std::abs(re), im), 1 - bits);
    const auto range = [&](const std::vector<std::size_t> &order, const std::vector<double> &key, double centre) {
        const auto first = std::lower_bound(order.begin(), order.end(), centre - window,
                                            [&](std::size_t k, double value) { return key[k] < value; });
        const auto last = std::upper_bound(first, order.end(), centre + window,
                                           [&](double value, std::size_t k) { return value < key[k]; });
        return std::pair(first, last);
    };
    const auto in_re = range(other.by_re, other.re, re);
    const auto in_im = range(other.by_im, other.im, im);
    const auto shorter = in_re.second - in_re.first <= in_im.second - in_im.first ? in_re : in_im;
    if (static_cast<std::size_t>(shorter.second - shorter.first) > scan_limit) {
        out.searched = false;
        return out;
    }
    for (auto it = shorter.first; it != shorter.second; ++it)
        consider(*it);
    return out;
}
} // namespace

std::vector<std::ptrdiff_t> conjugate_partners(std::span<const PlanePoint> points, int bits, std::size_t scan_limit) {
    std::vector<std::ptrdiff_t> partner(points.size(), -1);
    Side above, below;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const PlanePoint &p = points[i];
        if (!std::isfinite(p.re) || !std::isfinite(p.im) || p.im == 0)
            continue;
        (p.im > 0 ? above : below).add(i, p.re, std::abs(p.im));
    }
    if (above.size() == 0 || below.size() == 0)
        return partner;
    // Where every point of one side can be compared with every point of the other within scan_limit comparisons a
    // side, that is done: neither side then has more than scan_limit points, no window can hold more, the answer is
    // the sorted search's, and the sorts are saved (every low degree).
    const bool sorted = above.size() * below.size() > scan_limit;
    if (sorted) {
        above.sort();
        below.sort();
    }
    std::vector<Candidates> of_above(above.size()), of_below(below.size());
    for (std::size_t k = 0; k < above.size(); ++k)
        of_above[k] = candidates_of(above.re[k], above.im[k], below, sorted, bits, scan_limit);
    for (std::size_t k = 0; k < below.size(); ++k)
        of_below[k] = candidates_of(below.re[k], below.im[k], above, sorted, bits, scan_limit);
    // A pair: each searched, each with one candidate, the other. (The test is symmetric, so a point's one candidate
    // has that point among its own; where it has another too, or was not searched, there is no pair.)
    for (std::size_t k = 0; k < below.size(); ++k) {
        const Candidates &mine = of_below[k];
        if (!mine.searched || mine.count != 1)
            continue;
        const Candidates &theirs = of_above[mine.one];
        if (!theirs.searched || theirs.count != 1 || theirs.one != k)
            continue;
        partner[below.index[k]] = static_cast<std::ptrdiff_t>(above.index[mine.one]);
        partner[above.index[mine.one]] = static_cast<std::ptrdiff_t>(below.index[k]);
    }
    return partner;
}
} // namespace certpolysolve::detail
