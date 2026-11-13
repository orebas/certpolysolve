// FixedTransform (fixed_transform.hpp derives its headroom, bound and word size): the radix-2 discrete Fourier
// transform X_k = sum_j x_j w^(jk), w = exp(-2 pi i / N), N = 2^L, in fixed-point integers of K limbs, as balls that
// contain every X_k for every choice of the x_j in their balls, its error bound known before it runs (apply_into).
//
// apply_into and apply take the least K whose arithmetic bound is at most 2^(e - precision), 2^e above every midpoint
// component (limbs_for), or decline. FixedTransform::transform takes the roots exp(-2 pi i t / N), t < N / 2, from
// Impl::table (FLINT's unit roots held with G = 64 K - 2 fraction bits, built once a word size, unusable if a radius
// exceeds 2^-16 units). The conversion (to_fixed) holds each input component v as trunc(v 2^(F - e)),
// F = 64 K - L - headroom_bits, and sums the input radii into rho. The execution (dispatch, transform<K>) is
// decimation in frequency in place, each level taking A' = A + B and B' = floor((A - B) r 2^-G) for the pair's root
// r: the sum and difference exact (add, subtract), the product exact in 2 K limbs (rotate, multiply), the floor a
// shift (drop_root_scale). The enclosure (from_fixed) reads output k at its bit-reversed position (Impl::reversed),
// rounds it to the working precision (one ulp of radius if inexact) and adds arithmetic_bound (error_units_per_point
// N units) + rho in both coordinates. The file uses no floating point, so its results are the same on every CPU.
//
// State: Impl, the plan for one length: reversed, one root table a word size (tables, built on first use) and work.
// Entry points: the constructor and limbs_for from ModelTransform's constructor, apply_into from ModelTransform::apply
// (piecewise_models.cpp), under PiecewiseRootConfig::fixed_fft, which maintained_piecewise_config (piecewise_roots.hpp)
// sets. No other file in src/ calls apply, apply_with or arithmetic_bound.
// Not here: the choice of transform and FLINT's ball transform where this one declines (ModelTransform).
#include "numeric/transforms/fixed_transform.hpp"
#include "platform/parallel.hpp"
#include "core/telemetry.hpp"
#include "core/wide_int.hpp"
#include <algorithm>
#include <bit>
#include <span>
#include <stdexcept>
#include <vector>

namespace certpolysolve::detail {
namespace {
// FLINT's limb, so that its conversions take these arrays as they are.
using Limb = ulong;
using Wide = uint128;
constexpr int limb_bits = 64;
static_assert(FLINT_BITS == limb_bits, "the fixed transform is written for 64-bit limbs");

// Integers of K limbs, least significant first, in two's complement. A result may alias an operand.
template <int K>
inline void
add(Limb *r, const Limb *a, const Limb *b) {
    Wide carry = 0;
    for (int i = 0; i < K; ++i) {
        carry += static_cast<Wide>(a[i]) + b[i];
        r[i] = static_cast<Limb>(carry);
        carry >>= limb_bits;
    }
}
// a - (b where mask is all ones, zero where it is zero), as a + ~(b & mask) + 1.
template <int K>
inline void
subtract(Limb *r, const Limb *a, const Limb *b, Limb mask = ~Limb{0}) {
    Wide carry = 1;
    for (int i = 0; i < K; ++i) {
        carry += static_cast<Wide>(a[i]) + static_cast<Limb>(~(b[i] & mask));
        r[i] = static_cast<Limb>(carry);
        carry >>= limb_bits;
    }
}
// The exact product of two signed K-limb integers in 2 K limbs (r aliases neither). The unsigned product, then
// the two's complement correction: a signed operand is its unsigned reading minus 2^(64 K) when negative, so the
// other operand leaves the high half once for each negative one (modulo 2^(128 K)). No branch on the signs.
template <int K>
inline void
multiply(Limb *r, const Limb *a, const Limb *b) {
    for (int i = 0; i < 2 * K; ++i) r[i] = 0;
    for (int i = 0; i < K; ++i) {
        Limb carry = 0;
        for (int j = 0; j < K; ++j) {
            // At most (2^64 - 1)^2 + 2 (2^64 - 1) = 2^128 - 1.
            const Wide t = static_cast<Wide>(a[i]) * b[j] + r[i + j] + carry;
            r[i + j] = static_cast<Limb>(t);
            carry = static_cast<Limb>(t >> limb_bits);
        }
        r[i + K] = carry;
    }
    const Limb a_negative = Limb{0} - (a[K - 1] >> (limb_bits - 1));
    const Limb b_negative = Limb{0} - (b[K - 1] >> (limb_bits - 1));
    subtract<K>(r + K, r + K, b, a_negative);
    subtract<K>(r + K, r + K, a, b_negative);
}
// floor(p / 2^G) for a signed 2 K-limb p below 2^(2 G + 1) in size, G = 64 K - 2: its bits from G up, K limbs.
template <int K>
inline void
drop_root_scale(Limb *r, const Limb *p) {
    constexpr int spare = 2; // 64 K - G
    for (int i = 0; i < K; ++i) r[i] = (p[K - 1 + i] >> (limb_bits - spare)) | (p[K + i] << spare);
}
// (out_re + i out_im) = floor((d_re + i d_im) (c + i s) / 2^G) a component; the outputs alias no input. Always
// inlined into butterfly, as the transform's loop had it before W9 (below).
template <int K>
CERTPOLYSOLVE_WIDE_INT_ALWAYS_INLINE inline void
rotate(Limb *out_re, Limb *out_im, const Limb *d_re, const Limb *d_im, const Limb *c,
       const Limb *s) {
    Limb p[2 * K], q[2 * K];
    multiply<K>(p, d_re, c);
    multiply<K>(q, d_im, s);
    subtract<2 * K>(p, p, q);
    drop_root_scale<K>(out_re, p);
    multiply<K>(p, d_re, s);
    multiply<K>(q, d_im, c);
    add<2 * K>(p, p, q);
    drop_root_scale<K>(out_im, p);
}
// The transform of fixed_transform.hpp's header, in place on N = 2^L pairs of 2 K limbs (real part, then imaginary
// part); roots holds exp(-2 pi i t / 2^root_log2) for t below 2^(root_log2 - 1) in the same layout. The result is
// in bit-reversed order. The root of q = 0 is exactly one, and its product the difference itself.
// One butterfly: the pair at offsets q and q + half of the block at `start`, with its level's root stride. Always
// inlined: called from levels and from the split's chunks, Clang kept it out of line, and each butterfly paid a
// call (the W9 gate, 24 October 2026: 5.6% more instructions in the transform at one worker).
template <int K>
CERTPOLYSOLVE_WIDE_INT_ALWAYS_INLINE inline void
butterfly(Limb *z, const Limb *roots, size_t start, size_t half, size_t q, size_t stride) {
    constexpr size_t width = 2 * static_cast<size_t>(K);
    Limb *a = z + width * (start + q), *b = a + width * half;
    Limb d[2 * K];
    subtract<K>(d, a, b);
    subtract<K>(d + K, a + K, b + K);
    add<K>(a, a, b);
    add<K>(a + K, a + K, b + K);
    if (!q) {
        std::copy(d, d + 2 * K, b);
        return;
    }
    const Limb *r = roots + width * (q * stride);
    rotate<K>(b, b + K, d, d + K, r, r + K);
}
// Levels from <= level < to, on the blocks that lie in the points [low, high) (whole blocks at those levels).
template <int K>
void
levels(Limb *z, const Limb *roots, int root_log2, int L, int from, int to, size_t low, size_t high) {
    const size_t N = size_t{1} << L;
    for (int level = from; level < to; ++level) {
        const size_t block = N >> level, half = block / 2;
        const size_t stride = size_t{1} << (root_log2 - L + level);
        for (size_t start = low; start < high; start += block)
            for (size_t q = 0; q < half; ++q) butterfly<K>(z, roots, start, half, q, stride);
    }
}
// W9's split (parallel.hpp): a transform of at least 2^split_min_log2 points, on a solve with workers, takes its
// first L - sub_transform_log2 levels in chunks of butterfly_chunk butterflies, one parallel loop a level, and then
// the 2^(L - sub_transform_log2) sub-transforms of 2^sub_transform_log2 points those levels leave, one block each.
// Every butterfly is exact integer arithmetic on its own pair, and the pairs of a level are disjoint, so any split
// gives the serial transform's bits. The constants fix the blocks by the length alone: a level at 2,048 points takes
// about 27 microseconds (the W9 loop survey, L12), too little to share.
constexpr int split_min_log2 = 12, sub_transform_log2 = 11;
constexpr size_t butterfly_chunk = 1024;
template <int K>
void
transform(Limb *z, const Limb *roots, int root_log2, int L) {
    const size_t N = size_t{1} << L;
    if (L < split_min_log2 || parallel_width() < 2) {
        levels<K>(z, roots, root_log2, L, 0, L, 0, N);
        return;
    }
    const int split = L - sub_transform_log2;
    for (int level = 0; level < split; ++level) {
        const size_t block = N >> level, half = block / 2;
        const size_t stride = size_t{1} << (root_log2 - L + level);
        parallel_for(SplitSite::FixedTransform, N / 2 / butterfly_chunk, [&](size_t chunk) {
            for (size_t i = chunk * butterfly_chunk; i < (chunk + 1) * butterfly_chunk; ++i)
                butterfly<K>(z, roots, (i / half) * block, half, i % half, stride);
        });
    }
    const size_t M = size_t{1} << sub_transform_log2;
    parallel_for(SplitSite::FixedTransform, N / M,
                 [&](size_t sub) { levels<K>(z, roots, root_log2, L, split, L, sub * M, (sub + 1) * M); });
}
void
dispatch(int limbs, Limb *z, const Limb *roots, int root_log2, int L) {
    switch (limbs) {
    case 2: return transform<2>(z, roots, root_log2, L);
    case 3: return transform<3>(z, roots, root_log2, L);
    case 4: return transform<4>(z, roots, root_log2, L);
    case 5: return transform<5>(z, roots, root_log2, L);
    case 6: return transform<6>(z, roots, root_log2, L);
    case 7: return transform<7>(z, roots, root_log2, L);
    case 8: return transform<8>(z, roots, root_log2, L);
    default: throw std::logic_error("fixed transform word size out of range");
    }
}
int
ceil_log2(int x) {
    return std::bit_width(static_cast<unsigned>(x - 1));
}
// trunc(x 2^scale) as K limbs of two's complement, for a finite x with |x 2^scale| below 2^(64 K - 1). FLINT keeps
// a mantissa as normalised limbs (the top bit set) with x = 0.mantissa 2^exponent; _arf_get_integer_mpn writes the
// integer part of such a number, truncating its magnitude, in the ceil(exponent / 64) limbs it needs.
void
to_fixed(Limb *out, int K, const arf_t x, slong scale) {
    for (int i = 0; i < K; ++i) out[i] = 0;
    if (arf_is_zero(x)) return;
    nn_srcptr mantissa;
    slong limbs;
    ARF_GET_MPN_READONLY(mantissa, limbs, x);
    const slong exponent = fmpz_get_si(ARF_EXPREF(x)) + scale;
    if (exponent <= 0) return;
    _arf_get_integer_mpn(out, mantissa, limbs, exponent);
    if (ARF_SGNBIT(x)) {
        Wide carry = 1;
        for (int i = 0; i < K; ++i) {
            carry += static_cast<Limb>(~out[i]);
            out[i] = static_cast<Limb>(carry);
            carry >>= limb_bits;
        }
    }
}
// The K-limb two's complement integer `in` times 2^scale, as a ball of `precision` bits with the transform's
// bound: the magnitude through _arf_set_mpn_fixed, rounded toward zero as arb_set_round rounds, one unit in the
// last place added to the radius when that rounding was inexact.
void
from_fixed(arb_t x, const Limb *in, int K, slong scale, slong precision, const mag_t bound) {
    Limb magnitude[FixedTransform::maximum_limbs];
    const bool negative = (in[K - 1] >> (limb_bits - 1)) != 0;
    Wide carry = negative ? 1 : 0;
    for (int i = 0; i < K; ++i) {
        carry += negative ? static_cast<Limb>(~in[i]) : in[i];
        magnitude[i] = static_cast<Limb>(carry);
        carry >>= limb_bits;
    }
    const int inexact = _arf_set_mpn_fixed(arb_midref(x), magnitude, K, 0, negative, precision, ARB_RND);
    arf_mul_2exp_si(arb_midref(x), arb_midref(x), scale);
    if (inexact)
        arf_mag_set_ulp(arb_radref(x), arb_midref(x), precision);
    else
        mag_zero(arb_radref(x));
    mag_add(arb_radref(x), arb_radref(x), bound);
}
// The conversion's and the enclosure's loops over n points (each writes its own slot; apply_into's listed outputs are
// distinct): in fixed chunks of point_chunk points on a solve with workers, where n reaches 2^split_min_log2.
constexpr size_t point_chunk = 512;
template <class Body>
void
points(size_t n, const Body &body) {
    if (n < (size_t{1} << split_min_log2) || parallel_width() < 2) {
        for (size_t i = 0; i < n; ++i) body(i);
        return;
    }
    parallel_for(SplitSite::FixedTransform, (n + point_chunk - 1) / point_chunk, [&](size_t chunk) {
        for (size_t i = chunk * point_chunk; i < std::min(n, (chunk + 1) * point_chunk); ++i) body(i);
    });
}
// The bits a ball of a root of unity is computed with beyond its 64 K: its radius is then far below the 2^-16
// units the table accepts.
constexpr slong root_guard_bits = 64;
constexpr slong root_radius_units_log2 = -16;
} // namespace

int
FixedTransform::limbs_for(slong precision, int log2_length) {
    // Beyond the largest word's bits no word serves, and the sum below could pass slong: declined before anything
    // is added to the precision. The exponent's range is the constructor's.
    // History: docs/code-history-20261021.md, fixed_transform.cpp 1.
    if (precision < 1 || log2_length < 1 || precision > static_cast<slong>(limb_bits) * maximum_limbs) return 0;
    const slong bits = precision + 2 * static_cast<slong>(log2_length) + headroom_bits + ceil_log2(error_units_per_point);
    const slong limbs = std::max<slong>(minimum_limbs, (bits + limb_bits - 1) / limb_bits);
    return limbs <= maximum_limbs ? static_cast<int>(limbs) : 0;
}
void
FixedTransform::arithmetic_bound(mag_t bound, slong input_exponent, int log2_length, int limbs) {
    const slong fraction = static_cast<slong>(limb_bits) * limbs - log2_length - headroom_bits;
    mag_set_ui_2exp_si(bound, static_cast<ulong>(error_units_per_point), input_exponent - fraction + log2_length);
}

struct FixedTransform::Impl {
    int log2_length;
    slong length;
    // One table of roots a word size, built when first used.
    struct Table {
        std::vector<Limb> limbs;
        bool built = false, usable = false;
    };
    Table tables[maximum_limbs + 1];
    std::vector<Limb> work;
    Integer integer;
    // reversed[i]: the output held at position i of the transform's bit-reversed result.
    std::vector<uint32_t> reversed;

    // The exponent's range, checked before anything is shifted by it or sized from it (the constructor runs the
    // member initializer below on the checked value): a length of one needs no transform, and the positions of
    // a result are kept as 32-bit indices. History: docs/code-history-20261021.md, fixed_transform.cpp 2.
    static int checked(int exponent) {
        if (exponent < 1 || exponent > 31) throw std::invalid_argument("fixed transform length out of range");
        return exponent;
    }
    explicit Impl(int exponent)
      : log2_length(exponent)
      , length(slong{1} << exponent) {
        reversed.resize(static_cast<size_t>(length));
        for (size_t i = 0; i < reversed.size(); ++i) {
            uint32_t r = 0;
            for (int bit = 0; bit < exponent; ++bit) r |= static_cast<uint32_t>((i >> bit) & 1) << (exponent - 1 - bit);
            reversed[i] = r;
        }
    }

    const Table &table(int K) {
        Table &t = tables[K];
        if (t.built) return t;
        t.built = true;
        const slong G = static_cast<slong>(limb_bits) * K - 2;
        const slong count = std::max<slong>(1, length / 2);
        RootVector roots(count);
        _acb_vec_unit_roots(roots.data(), -length, count, static_cast<slong>(limb_bits) * K + root_guard_bits);
        t.limbs.assign(static_cast<size_t>(2 * K) * static_cast<size_t>(count), 0);
        for (slong i = 0; i < count; ++i)
            for (int part = 0; part < 2; ++part) {
                const arb_struct *x = part ? acb_imagref(roots[i]) : acb_realref(roots[i]);
                // The held component is the truncated midpoint: below one unit from it, and the midpoint within
                // the ball's radius of the root; twiddle_units = 2 covers both when the radius is small.
                if (!arb_is_finite(x) || mag_cmp_2exp_si(arb_radref(x), root_radius_units_log2 - G) > 0) return t;
                arf_get_fmpz_fixed_si(integer.value(), arb_midref(x), -G);
                fmpz_get_signed_ui_array(t.limbs.data() + static_cast<size_t>(2 * K) * static_cast<size_t>(i) +
                                             static_cast<size_t>(part * K),
                                         K, integer.value());
            }
        t.usable = true;
        return t;
    }
};
FixedTransform::FixedTransform(int log2_length)
  : impl_(std::make_unique<Impl>(Impl::checked(log2_length))) {}
FixedTransform::~FixedTransform() = default;

std::optional<EnclosedTransform>
FixedTransform::apply(const RootVector &input, slong precision, FixedTransformStatistics &stats) {
    const int limbs = limbs_for(precision, impl_->log2_length);
    if (!limbs) {
        ++stats.declined;
        return std::nullopt;
    }
    return apply_with(input, precision, limbs, stats);
}
std::optional<EnclosedTransform>
FixedTransform::apply_with(const RootVector &input, slong precision, int K, FixedTransformStatistics &stats) {
    if (K < minimum_limbs || K > maximum_limbs) throw std::invalid_argument("fixed transform word size out of range");
    EnclosedTransform enclosed(impl_->length);
    if (!transform(enclosed.values, input, precision, K, {}, stats)) return std::nullopt;
    return enclosed;
}
bool
FixedTransform::apply_into(RootVector &output, const RootVector &input, slong precision, std::span<const slong> wanted,
                           FixedTransformStatistics &stats) {
    const int limbs = limbs_for(precision, impl_->log2_length);
    if (!limbs) {
        ++stats.declined;
        return false;
    }
    return transform(output, input, precision, limbs, wanted, stats);
}
bool
FixedTransform::transform(RootVector &output, const RootVector &input, slong precision, int K,
                          std::span<const slong> wanted, FixedTransformStatistics &stats) {
    auto &p = *impl_;
    if (input.size() != p.length || output.size() != p.length)
        throw std::invalid_argument("fixed transform length does not match its plan");
    const auto declined = [&] {
        ++stats.declined;
        return false;
    };
    const int L = p.log2_length;
    const slong fraction = static_cast<slong>(limb_bits) * K - L - headroom_bits;
    if (precision < 2 || fraction < 1) return declined();
    const auto &table = p.table(K);
    if (!table.usable) return declined();
    const size_t width = 2 * static_cast<size_t>(K), N = static_cast<size_t>(p.length);
    slong exponent = 0;
    bool nonzero = false;
    Magnitude radii, bound;
    {
        StageTimer timer(stats.conversion_ns);
        // 2^exponent exceeds every midpoint component; the radii add up to rho.
        for (slong j = 0; j < p.length; ++j)
            for (int part = 0; part < 2; ++part) {
                const arb_struct *x = part ? acb_imagref(input[j]) : acb_realref(input[j]);
                if (!arb_is_finite(x)) return declined();
                mag_add(radii.value(), radii.value(), arb_radref(x));
                if (arf_is_zero(arb_midref(x))) continue;
                const slong size = arf_abs_bound_lt_2exp_si(arb_midref(x));
                // arf_abs_bound_lt_2exp_si clamps at ARF_PREC_EXACT: an exponent near it is not a real size.
                if (size >= ARF_PREC_EXACT / 2 || size <= -ARF_PREC_EXACT / 2) return declined();
                exponent = nonzero ? std::max(exponent, size) : size;
                nonzero = true;
            }
        p.work.resize(width * N);
        // trunc(v 2^(F - e)): below 2^F in size, within one unit of v. Each point writes its own slot.
        const auto convert = [&](size_t j) {
            Limb *slot = p.work.data() + width * j;
            to_fixed(slot, K, arb_midref(acb_realref(input[static_cast<slong>(j)])), fraction - exponent);
            to_fixed(slot + K, K, arb_midref(acb_imagref(input[static_cast<slong>(j)])), fraction - exponent);
        };
        points(N, convert);
    }
    {
        StageTimer timer(stats.execution_ns);
        dispatch(K, p.work.data(), table.limbs.data(), L, L);
    }
    StageTimer timer(stats.enclosure_ns);
    arithmetic_bound(bound.value(), exponent, L, K);
    mag_add(bound.value(), bound.value(), radii.value());
    if (!mag_is_finite(bound.value())) return declined();
    // Position i of the bit-reversed result holds output reversed[i]; reversal is its own inverse, so output k
    // is at position reversed[k]. Rounding a midpoint to the working precision adds its own error to the radius.
    const auto publish = [&](slong k) {
        const Limb *slot = p.work.data() + width * static_cast<size_t>(p.reversed[static_cast<size_t>(k)]);
        from_fixed(acb_realref(output[k]), slot, K, exponent - fraction, precision, bound.value());
        from_fixed(acb_imagref(output[k]), slot + K, K, exponent - fraction, precision, bound.value());
    };
    if (wanted.empty()) {
        points(N, [&](size_t k) { publish(static_cast<slong>(k)); });
    } else {
        for (const slong k : wanted)
            if (k < 0 || k >= p.length) throw std::out_of_range("fixed transform output index");
        points(wanted.size(), [&](size_t i) { publish(wanted[i]); });
    }
    ++stats.accepted;
    stats.limbs = static_cast<uint64_t>(K);
    return true;
}
} // namespace certpolysolve::detail
