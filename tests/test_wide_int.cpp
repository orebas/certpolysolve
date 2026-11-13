// The portable 128-bit integers and checked helpers of src/core/wide_int.hpp (the library's arithmetic where the compiler
// has no __int128 and no overflow built-ins, MSVC) against the native types and built-ins where they exist (GCC,
// Clang), and everywhere against a model in 32-bit words: every operator and conversion the sites use, on edge
// values (0, 1, -1, the powers of two from 2^0 to 2^127 and their neighbours, every shift amount 0..127, operands
// of every sign) and on random operands from a fixed seed.
#include "core/wide_int.hpp"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <type_traits>
#include <vector>

namespace d = certpolysolve::detail;
using d::portable::Int128;
using d::portable::UInt128;
#ifdef __SIZEOF_INT128__
__extension__ using NativeU = unsigned __int128;
__extension__ using NativeS = __int128;
#endif

static_assert(d::uint128_bits == 128 && d::wide_limb_bits == 64, "two limbs of 64 bits");
static_assert(sizeof(UInt128) == 16 && sizeof(Int128) == 16, "the portable types are two limbs");
#ifdef CERTPOLYSOLVE_PORTABLE_INT128
static_assert(std::is_same_v<d::uint128, UInt128> && std::is_same_v<d::int128, Int128>,
              "CERTPOLYSOLVE_PORTABLE_INT128 selects the portable classes");
#elif defined(__SIZEOF_INT128__)
static_assert(std::is_same_v<d::uint128, NativeU> && std::is_same_v<d::int128, NativeS>,
              "a compiler with __int128 uses it unless the fallback is forced");
#else
static_assert(std::is_same_v<d::uint128, UInt128> && std::is_same_v<d::int128, Int128>,
              "a compiler without __int128 (MSVC) uses the portable classes");
#endif
// Constant evaluation, as the sites' constants (Signed{1} << (product_bits - 1)) need it.
static_assert((UInt128(1) << 127) == UInt128::from_limbs(0, std::uint64_t{1} << 63), "a shift into the high limb");
static_assert((Int128(-1) >> 100) == Int128(-1), "an arithmetic shift keeps the sign");
static_assert((Int128(1) << 125) > Int128(0) && (Int128(1) << 127) < Int128(0), "bit 127 is the sign");
static_assert(static_cast<std::uint64_t>(UInt128::from_limbs(5, 7)) == 5, "a conversion keeps the low limb");
static_assert(UInt128(-1) == UInt128::from_limbs(~std::uint64_t{0}, ~std::uint64_t{0}), "sign extension");
static_assert(UInt128(std::uint64_t{1} << 63) * UInt128(4) == UInt128::from_limbs(0, 2), "a product modulo 2^128");

namespace {
int failures = 0;
void check(bool ok, const char *what, std::uint64_t low = 0, std::uint64_t high = 0, int k = -1) {
    if (ok) return;
    ++failures;
    if (failures <= 40) std::fprintf(stderr, "FAIL %s (operand low %#llx high %#llx, k %d)\n", what,
                                     static_cast<unsigned long long>(low), static_cast<unsigned long long>(high), k);
}

// The model: four 32-bit words, least significant first, two's complement.
struct Words {
    std::uint32_t w[4] = {0, 0, 0, 0};
    static Words from(std::uint64_t low, std::uint64_t high) {
        Words r;
        r.w[0] = static_cast<std::uint32_t>(low);
        r.w[1] = static_cast<std::uint32_t>(low >> 32);
        r.w[2] = static_cast<std::uint32_t>(high);
        r.w[3] = static_cast<std::uint32_t>(high >> 32);
        return r;
    }
    std::uint64_t low() const { return w[0] | (std::uint64_t{w[1]} << 32); }
    std::uint64_t high() const { return w[2] | (std::uint64_t{w[3]} << 32); }
    bool bit(int i) const { return (w[i / 32] >> (i % 32)) & 1; }
    bool sign() const { return bit(127); }
    bool operator==(const Words &o) const { return w[0] == o.w[0] && w[1] == o.w[1] && w[2] == o.w[2] && w[3] == o.w[3]; }
    static Words add(Words a, Words b) {
        Words r;
        std::uint64_t carry = 0;
        for (int i = 0; i < 4; ++i) {
            carry += std::uint64_t{a.w[i]} + b.w[i];
            r.w[i] = static_cast<std::uint32_t>(carry);
            carry >>= 32;
        }
        return r;
    }
    static Words bit_not(Words a) {
        for (auto &x : a.w) x = ~x;
        return a;
    }
    static Words neg(Words a) { return add(bit_not(a), from(1, 0)); }
    static Words sub(Words a, Words b) { return add(a, neg(b)); }
    static Words mul(Words a, Words b) {
        std::uint32_t r[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        for (int i = 0; i < 4; ++i) {
            std::uint64_t carry = 0;
            for (int j = 0; j < 4; ++j) {
                carry += std::uint64_t{a.w[i]} * b.w[j] + r[i + j];
                r[i + j] = static_cast<std::uint32_t>(carry);
                carry >>= 32;
            }
            r[i + 4] = static_cast<std::uint32_t>(carry);
        }
        Words out;
        for (int i = 0; i < 4; ++i) out.w[i] = r[i];
        return out;
    }
    static Words bit_and(Words a, Words b) {
        for (int i = 0; i < 4; ++i) a.w[i] &= b.w[i];
        return a;
    }
    static Words bit_or(Words a, Words b) {
        for (int i = 0; i < 4; ++i) a.w[i] |= b.w[i];
        return a;
    }
    static Words bit_xor(Words a, Words b) {
        for (int i = 0; i < 4; ++i) a.w[i] ^= b.w[i];
        return a;
    }
    // Bit by bit: slow and obviously right.
    static Words shl(Words a, int k) {
        Words r;
        for (int i = 0; i < 128; ++i)
            if (i - k >= 0 && a.bit(i - k)) r.w[i / 32] |= std::uint32_t{1} << (i % 32);
        return r;
    }
    static Words shr(Words a, int k, bool arithmetic) {
        Words r;
        for (int i = 0; i < 128; ++i) {
            const bool bit = i + k < 128 ? a.bit(i + k) : arithmetic && a.sign();
            if (bit) r.w[i / 32] |= std::uint32_t{1} << (i % 32);
        }
        return r;
    }
    static int compare_unsigned(Words a, Words b) {
        for (int i = 3; i >= 0; --i)
            if (a.w[i] != b.w[i]) return a.w[i] < b.w[i] ? -1 : 1;
        return 0;
    }
    static int compare_signed(Words a, Words b) {
        if (a.sign() != b.sign()) return a.sign() ? -1 : 1;
        return compare_unsigned(a, b);
    }
    static Words from_signed(std::int64_t v) { return from(static_cast<std::uint64_t>(v), v < 0 ? ~std::uint64_t{0} : 0); }
};
Words words(UInt128 x) { return Words::from(x.low(), x.high()); }
Words words(Int128 x) { return Words::from(x.low(), x.high()); }
bool same(UInt128 x, Words w) { return words(x) == w; }
bool same(Int128 x, Words w) { return words(x) == w; }

#ifdef __SIZEOF_INT128__
NativeU native_unsigned(std::uint64_t low, std::uint64_t high) { return (static_cast<NativeU>(high) << 64) | low; }
NativeS native_signed(std::uint64_t low, std::uint64_t high) { return static_cast<NativeS>(native_unsigned(low, high)); }
bool same(UInt128 x, NativeU n) { return x.low() == static_cast<std::uint64_t>(n) && x.high() == static_cast<std::uint64_t>(n >> 64); }
bool same(Int128 x, NativeS n) { return same(UInt128(x), static_cast<NativeU>(n)); }
bool same(Words w, NativeU n) { return w.low() == static_cast<std::uint64_t>(n) && w.high() == static_cast<std::uint64_t>(n >> 64); }
#endif

// One operand's bit pattern.
struct Operand {
    std::uint64_t low, high;
};
std::vector<Operand> edge_operands() {
    std::vector<Operand> out;
    const auto add = [&](std::uint64_t low, std::uint64_t high) { out.push_back({low, high}); };
    add(0, 0);
    add(1, 0);
    add(~std::uint64_t{0}, ~std::uint64_t{0}); // -1
    add(2, 0);
    add(3, 0);
    add(0x5555555555555555, 0x5555555555555555);
    add(0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA);
    add(0x0123456789ABCDEF, 0xFEDCBA9876543210);
    // 2^k, 2^k - 1, 2^k + 1 and -2^k for every k: 2^63, 2^64 - 1, 2^64, 2^127 and -2^127 among them.
    for (int k = 0; k < 128; ++k) {
        const Words power = Words::shl(Words::from(1, 0), k);
        for (const Words v : {power, Words::sub(power, Words::from(1, 0)), Words::add(power, Words::from(1, 0)), Words::neg(power)})
            add(v.low(), v.high());
    }
    return out;
}
std::vector<Operand> random_operands(std::mt19937_64 &random, size_t count) {
    std::vector<Operand> out;
    out.reserve(count);
    while (out.size() < count) {
        std::uint64_t low = random(), high = random();
        switch (random() % 8) {
        case 0: break;                                                         // uniform
        case 1: high = 0; break;                                               // one limb
        case 2: low = 0; break;                                                // the high limb alone
        case 3: {                                                              // a small signed value
            const auto v = static_cast<std::int64_t>(random() % 2001) - 1000;
            low = static_cast<std::uint64_t>(v);
            high = v < 0 ? ~std::uint64_t{0} : 0;
            break;
        }
        case 4: {                                                              // 2^k plus a small signed offset
            const Words power = Words::shl(Words::from(1, 0), static_cast<int>(random() % 128));
            const Words v = Words::add(power, Words::from_signed(static_cast<std::int64_t>(random() % 7) - 3));
            low = v.low();
            high = v.high();
            break;
        }
        case 5: {                                                              // a random width
            const Words v = Words::shr(Words::from(low, high), static_cast<int>(random() % 128), false);
            low = v.low();
            high = v.high();
            break;
        }
        case 6: high = (low >> 63) ? ~std::uint64_t{0} : 0; break;              // a sign-extended 64-bit value
        default: {                                                             // the negation of a uniform value
            const Words v = Words::neg(Words::from(low, high));
            low = v.low();
            high = v.high();
            break;
        }
        }
        out.push_back({low, high});
    }
    return out;
}

template <class T>
bool bits_equal(T a, T b) {
    return a == b && std::signbit(a) == std::signbit(b);
}

void unary_cases(Operand o, const std::vector<int> &shifts) {
    const UInt128 u = UInt128::from_limbs(o.low, o.high);
    const Int128 s = Int128::from_limbs(o.low, o.high);
    const Words w = Words::from(o.low, o.high);
    check(same(-u, Words::neg(w)), "unsigned negation", o.low, o.high);
    check(same(~u, Words::bit_not(w)), "unsigned complement", o.low, o.high);
    check(same(-s, Words::neg(w)), "signed negation", o.low, o.high);
    check(same(~s, Words::bit_not(w)), "signed complement", o.low, o.high);
    check(static_cast<bool>(u) == (o.low != 0 || o.high != 0), "unsigned bool", o.low, o.high);
    check(static_cast<bool>(s) == (o.low != 0 || o.high != 0), "signed bool", o.low, o.high);
    check(!u == (o.low == 0 && o.high == 0), "unsigned not", o.low, o.high);
    check(static_cast<std::uint64_t>(u) == o.low, "to uint64", o.low, o.high);
    check(static_cast<std::int64_t>(s) == static_cast<std::int64_t>(o.low), "to int64", o.low, o.high);
    check(static_cast<std::uint32_t>(u) == static_cast<std::uint32_t>(o.low), "to uint32", o.low, o.high);
    check(static_cast<int>(s) == static_cast<int>(static_cast<std::uint32_t>(o.low)), "to int", o.low, o.high);
    check(static_cast<unsigned long>(u) == static_cast<unsigned long>(o.low), "to unsigned long", o.low, o.high);
    check(static_cast<long>(s) == static_cast<long>(o.low), "to long", o.low, o.high);
    // The casts between the two readings keep the bits.
    check(UInt128(s).low() == o.low && UInt128(s).high() == o.high, "signed to unsigned", o.low, o.high);
    check(Int128(u).low() == o.low && Int128(u).high() == o.high, "unsigned to signed", o.low, o.high);
    check(s.negative() == w.sign(), "sign", o.low, o.high);
    check((s < 0) == w.sign(), "sign by comparison", o.low, o.high);
    for (const int k : shifts) {
        check(same(u << k, Words::shl(w, k)), "unsigned <<", o.low, o.high, k);
        check(same(u >> k, Words::shr(w, k, false)), "unsigned >>", o.low, o.high, k);
        check(same(s << k, Words::shl(w, k)), "signed <<", o.low, o.high, k);
        check(same(s >> k, Words::shr(w, k, true)), "signed >>", o.low, o.high, k);
        // Amounts as long, the type of the sites' shift variables.
        check(same(u >> static_cast<long>(k), Words::shr(w, k, false)), "unsigned >> long", o.low, o.high, k);
        check(same(s >> static_cast<long>(k), Words::shr(w, k, true)), "signed >> long", o.low, o.high, k);
        UInt128 cu = u;
        cu <<= k;
        check(same(cu, Words::shl(w, k)), "unsigned <<=", o.low, o.high, k);
        cu = u;
        cu >>= k;
        check(same(cu, Words::shr(w, k, false)), "unsigned >>=", o.low, o.high, k);
        Int128 cs = s;
        cs >>= k;
        check(same(cs, Words::shr(w, k, true)), "signed >>=", o.low, o.high, k);
    }
    // Shifts beyond the native range move every bit out (the native types leave them undefined).
    check((u << 128) == UInt128(0) && (u >> 128) == UInt128(0), "unsigned shift by 128", o.low, o.high);
    check((s >> 128) == (w.sign() ? Int128(-1) : Int128(0)) && (s << 200) == Int128(0), "signed shift by 128+", o.low, o.high);
    // Floating conversions: exact below 2^53 (double) and 2^64 (an 80-bit long double) by construction.
    if (o.high == 0 && o.low < (std::uint64_t{1} << 53)) {
        check(static_cast<double>(u) == static_cast<double>(o.low), "to double, exact", o.low, o.high);
        check(static_cast<long double>(u) == static_cast<long double>(o.low), "to long double, exact", o.low, o.high);
    }
    if (o.high == ~std::uint64_t{0} && o.low >= ~(std::uint64_t{1} << 53))
        check(static_cast<double>(s) == static_cast<double>(static_cast<std::int64_t>(o.low)), "negative to double, exact", o.low, o.high);
#ifdef __SIZEOF_INT128__
    const NativeU nu = native_unsigned(o.low, o.high);
    const NativeS ns = native_signed(o.low, o.high);
    check(same(w, nu), "the model's bits", o.low, o.high);
    check(same(-u, static_cast<NativeU>(-nu)), "unsigned negation, native", o.low, o.high);
    check(same(~u, ~nu), "unsigned complement, native", o.low, o.high);
    check(same(-s, static_cast<NativeS>(-static_cast<NativeU>(ns))), "signed negation, native", o.low, o.high);
    check(static_cast<bool>(u) == static_cast<bool>(nu), "unsigned bool, native", o.low, o.high);
    check(static_cast<std::int64_t>(s) == static_cast<std::int64_t>(ns), "to int64, native", o.low, o.high);
    check(static_cast<std::uint64_t>(s) == static_cast<std::uint64_t>(ns), "signed to uint64, native", o.low, o.high);
    check(static_cast<int>(s) == static_cast<int>(ns), "to int, native", o.low, o.high);
    check(same(UInt128(s), static_cast<NativeU>(ns)), "signed to unsigned, native", o.low, o.high);
    check((s < 0) == (ns < 0), "sign, native", o.low, o.high);
    for (const int k : shifts) {
        check(same(u << k, nu << k), "unsigned <<, native", o.low, o.high, k);
        check(same(u >> k, nu >> k), "unsigned >>, native", o.low, o.high, k);
        check(same(s >> k, ns >> k), "signed >>, native", o.low, o.high, k);
        // A left shift of a negative or overflowing signed value is defined on the unsigned reading (C++20).
        check(same(s << k, static_cast<NativeS>(static_cast<NativeU>(ns) << k)), "signed <<, native", o.low, o.high, k);
    }
    check(bits_equal(static_cast<long double>(u), static_cast<long double>(nu)), "unsigned to long double, native", o.low, o.high);
    check(bits_equal(static_cast<double>(u), static_cast<double>(nu)), "unsigned to double, native", o.low, o.high);
    check(bits_equal(static_cast<float>(u), static_cast<float>(nu)), "unsigned to float, native", o.low, o.high);
    check(bits_equal(static_cast<long double>(s), static_cast<long double>(ns)), "signed to long double, native", o.low, o.high);
    check(bits_equal(static_cast<double>(s), static_cast<double>(ns)), "signed to double, native", o.low, o.high);
    check(bits_equal(static_cast<float>(s), static_cast<float>(ns)), "signed to float, native", o.low, o.high);
#endif
}

void binary_cases(Operand a, Operand b) {
    const UInt128 ua = UInt128::from_limbs(a.low, a.high), ub = UInt128::from_limbs(b.low, b.high);
    const Int128 sa = Int128::from_limbs(a.low, a.high), sb = Int128::from_limbs(b.low, b.high);
    const Words wa = Words::from(a.low, a.high), wb = Words::from(b.low, b.high);
    check(same(ua + ub, Words::add(wa, wb)), "unsigned +", a.low, a.high);
    check(same(ua - ub, Words::sub(wa, wb)), "unsigned -", a.low, a.high);
    check(same(ua * ub, Words::mul(wa, wb)), "unsigned *", a.low, a.high);
    check(same(ua & ub, Words::bit_and(wa, wb)), "unsigned &", a.low, a.high);
    check(same(ua | ub, Words::bit_or(wa, wb)), "unsigned |", a.low, a.high);
    check(same(ua ^ ub, Words::bit_xor(wa, wb)), "unsigned ^", a.low, a.high);
    check(same(sa + sb, Words::add(wa, wb)), "signed +", a.low, a.high);
    check(same(sa - sb, Words::sub(wa, wb)), "signed -", a.low, a.high);
    check(same(sa * sb, Words::mul(wa, wb)), "signed *", a.low, a.high);
    check(same(sa & sb, Words::bit_and(wa, wb)), "signed &", a.low, a.high);
    check(same(sa | sb, Words::bit_or(wa, wb)), "signed |", a.low, a.high);
    check(same(sa ^ sb, Words::bit_xor(wa, wb)), "signed ^", a.low, a.high);
    const int cu = Words::compare_unsigned(wa, wb), cs = Words::compare_signed(wa, wb);
    check((ua == ub) == (cu == 0) && (ua != ub) == (cu != 0), "unsigned ==", a.low, a.high);
    check((ua < ub) == (cu < 0) && (ua <= ub) == (cu <= 0) && (ua > ub) == (cu > 0) && (ua >= ub) == (cu >= 0), "unsigned <", a.low, a.high);
    check((sa == sb) == (cs == 0) && (sa != sb) == (cs != 0), "signed ==", a.low, a.high);
    check((sa < sb) == (cs < 0) && (sa <= sb) == (cs <= 0) && (sa > sb) == (cs > 0) && (sa >= sb) == (cs >= 0), "signed <", a.low, a.high);
    // Compound assignments, and operands of 64 bits converted at the operator (the sites' Wide + limb).
    UInt128 c = ua;
    c += ub;
    check(same(c, Words::add(wa, wb)), "unsigned +=", a.low, a.high);
    c = ua;
    c -= ub;
    check(same(c, Words::sub(wa, wb)), "unsigned -=", a.low, a.high);
    c = ua;
    c *= ub;
    check(same(c, Words::mul(wa, wb)), "unsigned *=", a.low, a.high);
    c = ua;
    c |= ub;
    check(same(c, Words::bit_or(wa, wb)), "unsigned |=", a.low, a.high);
    check(same(ua + b.low, Words::add(wa, Words::from(b.low, 0))), "unsigned + uint64", a.low, a.high);
    check(same(ua * b.low, Words::mul(wa, Words::from(b.low, 0))), "unsigned * uint64", a.low, a.high);
    check(same(ua | b.low, Words::bit_or(wa, Words::from(b.low, 0))), "unsigned | uint64", a.low, a.high);
    const auto signed_low = static_cast<std::int64_t>(b.low);
    check(same(sa - signed_low, Words::sub(wa, Words::from_signed(signed_low))), "signed - int64", a.low, a.high);
    check(same(sa * signed_low, Words::mul(wa, Words::from_signed(signed_low))), "signed * int64", a.low, a.high);
    Int128 d = sa;
    d -= sb;
    check(same(d, Words::sub(wa, wb)), "signed -=", a.low, a.high);
    // Division and remainder (the divisor not zero, as for the native type): a = q b + r with r < b, and a divisor of
    // 64 bits converted at the operator (the sites' Wide / limb).
    if (ub != UInt128(0)) {
        const UInt128 q = ua / ub, r = ua % ub;
        check(r < ub && same(q * ub + r, wa), "unsigned / and %", a.low, a.high);
        UInt128 e = ua;
        e /= ub;
        check(e == q, "unsigned /=", a.low, a.high);
        e = ua;
        e %= ub;
        check(e == r, "unsigned %=", a.low, a.high);
    }
    if (b.low != 0)
        check(ua / b.low == ua / UInt128(b.low) && ua % b.low == ua % UInt128(b.low), "unsigned / uint64", a.low,
              a.high);
#ifdef __SIZEOF_INT128__
    const NativeU na = native_unsigned(a.low, a.high), nb = native_unsigned(b.low, b.high);
    const NativeS ma = native_signed(a.low, a.high), mb = native_signed(b.low, b.high);
    check(same(ua + ub, static_cast<NativeU>(na + nb)), "unsigned +, native", a.low, a.high);
    check(same(ua - ub, static_cast<NativeU>(na - nb)), "unsigned -, native", a.low, a.high);
    check(same(ua * ub, static_cast<NativeU>(na * nb)), "unsigned *, native", a.low, a.high);
    check((ua < ub) == (na < nb) && (ua == ub) == (na == nb) && (ua >= ub) == (na >= nb), "unsigned <, native", a.low, a.high);
    if (nb != 0)
        check(same(ua / ub, static_cast<NativeU>(na / nb)) && same(ua % ub, static_cast<NativeU>(na % nb)),
              "unsigned / and %, native", a.low, a.high);
    check((sa < sb) == (ma < mb) && (sa == sb) == (ma == mb) && (sa >= sb) == (ma >= mb), "signed <, native", a.low, a.high);
    // Signed sums and products on the unsigned reading: the native signed operation is undefined on overflow.
    check(same(sa + sb, static_cast<NativeS>(static_cast<NativeU>(ma) + static_cast<NativeU>(mb))), "signed +, native", a.low, a.high);
    check(same(sa - sb, static_cast<NativeS>(static_cast<NativeU>(ma) - static_cast<NativeU>(mb))), "signed -, native", a.low, a.high);
    check(same(sa * sb, static_cast<NativeS>(static_cast<NativeU>(ma) * static_cast<NativeU>(mb))), "signed *, native", a.low, a.high);
    // Where the product fits, the signed native product itself.
    if (a.high == 0 && a.low < (std::uint64_t{1} << 62) && b.high == 0 && b.low < (std::uint64_t{1} << 62))
        check(same(sa * sb, ma * mb), "signed * in range, native", a.low, a.high);
#endif
}

// Constructions from the integer types, with sign extension, against the model.
void construction_cases(std::uint64_t bits) {
    const auto s64 = static_cast<std::int64_t>(bits);
    const auto s32 = static_cast<std::int32_t>(bits);
    const auto u32 = static_cast<std::uint32_t>(bits);
    check(same(UInt128(bits), Words::from(bits, 0)), "unsigned from uint64", bits);
    check(same(UInt128(s64), Words::from_signed(s64)), "unsigned from int64", bits);
    check(same(UInt128(s32), Words::from_signed(s32)), "unsigned from int32", bits);
    check(same(UInt128(u32), Words::from(u32, 0)), "unsigned from uint32", bits);
    check(same(Int128(bits), Words::from(bits, 0)), "signed from uint64", bits);
    check(same(Int128(s64), Words::from_signed(s64)), "signed from int64", bits);
    check(same(Int128(s32), Words::from_signed(s32)), "signed from int32", bits);
    check(same(Int128(static_cast<long>(s64)), Words::from_signed(static_cast<long>(s64))), "signed from long", bits);
    check(same(Int128(static_cast<unsigned long>(bits)), Words::from(static_cast<unsigned long>(bits), 0)), "signed from unsigned long", bits);
    UInt128 assigned = 0;
    assigned = bits;
    check(same(assigned, Words::from(bits, 0)), "unsigned assigned from uint64", bits);
    Int128 initialised = s64;
    check(same(initialised, Words::from_signed(s64)), "signed copy-initialised from int64", bits);
    check(UInt128(s64) == UInt128(Int128(s64)) && (Int128(s64) == s64) && (UInt128(bits) == bits), "comparisons with integers", bits);
#ifdef __SIZEOF_INT128__
    check(same(UInt128(s64), static_cast<NativeU>(s64)), "unsigned from int64, native", bits);
    check(same(Int128(s32), static_cast<NativeS>(s32)), "signed from int32, native", bits);
#endif
}

// The checked helpers: the selected ones (the built-ins on GCC and Clang unless the fallback is forced) and
// the portable ones, against the model's exact product and difference, and against the built-ins.
template <class T>
void checked_cases(T a, T b) {
    const Words exact_product = Words::mul(Words::from_signed(a), Words::from_signed(b));
    const Words exact_difference = Words::sub(Words::from_signed(a), Words::from_signed(b));
    const auto fits = [](Words w) {
        const Words low = Words::from_signed(static_cast<T>(static_cast<std::int64_t>(w.low())));
        return w == low; // the value is the sign extension of its low bits in T
    };
    const auto wrapped = [](Words w) { return static_cast<T>(static_cast<std::int64_t>(w.low())); };
    T p = 0, q = 0, r = 0, s = 0;
    const bool mo = d::multiply_overflows(a, b, &p), po = d::portable::multiply_overflows(a, b, &q);
    const bool so = d::subtract_overflows(a, b, &r), qo = d::portable::subtract_overflows(a, b, &s);
    check(mo == !fits(exact_product) && p == wrapped(exact_product), "multiply_overflows", static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(b));
    check(po == !fits(exact_product) && q == wrapped(exact_product), "portable multiply_overflows", static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(b));
    check(so == !fits(exact_difference) && r == wrapped(exact_difference), "subtract_overflows", static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(b));
    check(qo == !fits(exact_difference) && s == wrapped(exact_difference), "portable subtract_overflows", static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(b));
    // The output aliasing an input, as machine_adapter.hpp writes it.
    T aliased = a;
    const bool ao = d::subtract_overflows(aliased, b, &aliased);
    check(ao == so && aliased == r, "subtract_overflows with an aliased output", static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(b));
#if defined(__has_builtin)
#if __has_builtin(__builtin_mul_overflow) && __has_builtin(__builtin_sub_overflow)
    T bp = 0, bs = 0;
    const bool bmo = __builtin_mul_overflow(a, b, &bp), bso = __builtin_sub_overflow(a, b, &bs);
    check(bmo == po && bp == q, "portable multiply_overflows against the built-in", static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(b));
    check(bso == qo && bs == s, "portable subtract_overflows against the built-in", static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(b));
#endif
#endif
}
template <class T>
void checked_suite(std::mt19937_64 &random) {
    constexpr T minimum = std::numeric_limits<T>::min(), maximum = std::numeric_limits<T>::max();
    std::vector<T> edges = {0, 1, -1, 2, -2, 3, -3, minimum, maximum, static_cast<T>(minimum + 1), static_cast<T>(maximum - 1),
                            static_cast<T>(minimum / 2), static_cast<T>(maximum / 2), static_cast<T>(minimum / 2 - 1), static_cast<T>(maximum / 2 + 1)};
    // The powers of two of the type and their neighbours: the boundaries of a product that fits.
    for (int k = 0; k < std::numeric_limits<T>::digits; ++k) {
        const T power = static_cast<T>(T{1} << k);
        for (const T v : {power, static_cast<T>(power - 1), static_cast<T>(power + 1)})
            for (const T w : {v, static_cast<T>(-v)}) edges.push_back(w);
    }
    for (const T a : edges)
        for (const T b : edges) checked_cases(a, b);
    for (int i = 0; i < 100000; ++i) {
        const auto pick = [&] {
            const std::uint64_t bits = random();
            switch (random() % 4) {
            case 0: return static_cast<T>(bits);
            case 1: return static_cast<T>(static_cast<std::int64_t>(bits) >> (random() % 64));
            case 2: return static_cast<T>(static_cast<std::int64_t>(bits % 1000) - 500);
            default: return edges[bits % edges.size()];
            }
        };
        checked_cases(pick(), pick());
    }
}
} // namespace

int main() {
    std::mt19937_64 random(20261020);
    std::vector<int> every_shift;
    for (int k = 0; k < 128; ++k) every_shift.push_back(k);
    const auto edges = edge_operands();
    for (const Operand o : edges) unary_cases(o, every_shift);
    for (const Operand a : edges)
        for (const Operand b : edges) binary_cases(a, b);
    const auto randoms = random_operands(random, 60000);
    for (const Operand o : randoms) {
        const std::vector<int> shifts = {static_cast<int>(random() % 128), 0, 1, 63, 64, 65, 127};
        unary_cases(o, shifts);
    }
    for (size_t i = 0; i + 1 < randoms.size(); ++i) binary_cases(randoms[i], randoms[i + 1]);
    for (size_t i = 0; i < randoms.size(); i += 7) binary_cases(randoms[i], edges[i % edges.size()]);
    for (const Operand o : edges) construction_cases(o.low);
    for (int i = 0; i < 20000; ++i) construction_cases(random());
    checked_suite<std::int64_t>(random);
    checked_suite<long>(random);
    checked_suite<long long>(random);
    checked_suite<int>(random);
    if (failures) {
        std::fprintf(stderr, "test_wide_int: %d failures\n", failures);
        return 1;
    }
    std::printf("test_wide_int: %zu edge operands, %zu random operands, every operator and conversion of the sites agree\n",
                edges.size(), randoms.size());
    return 0;
}
