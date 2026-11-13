#pragma once
#include "core/raii.hpp"
#include "core/wide_int.hpp"
#include <cstdint>
#include <span>
#include <vector>

namespace certpolysolve::detail {
// The node products of a secular model, prod_{j != i} (b_i - b_j), for nodes whose parts are short numbers (a
// binary64 number, or at most two limbs), in integer mantissas with a relative error bound that is known before
// a product is formed.
//
// Why: a regeneration needs one such product a node, n (n - 1) differences and multiplications in all. In every
// regeneration but the last the nodes are binary64 numbers (they come from the binary64 packets). In ball
// arithmetic at 117 bits a factor is a ball subtraction, a magnitude and a ball multiplication, and the products
// are a large part of a hard secular solve; the integers below take a fraction of that time. The nodes of a last
// regeneration, and the ones a certificate attempt's Newton steps hand over, have up to 128 bits: their factors
// are formed from 256-bit differences, and the product as for binary64 nodes.
// History: docs/code-history-20261021.md, node_products.hpp 1.
//
// A factor between binary64 parts. A part is m 2^q with |m| < 2^53 an integer. The difference of two parts is
// formed as an integer D with an exponent q0 = max(q_small, q_large - difference_shift): exact when the
// exponents are within difference_shift of each other, and otherwise the part of smaller exponent loses less
// than one unit 2^q0 to a shift, while |D| is then at least 2^(52 + difference_shift - 1). The two parts'
// differences are brought to one exponent f by which the larger has its leading bit at position factor_bits (a
// further truncation below one unit each), so a held factor d^ has its larger part in
// [2^(factor_bits - 1), 2^factor_bits) and |d^ - d| < 2 sqrt(2) units: relative error at most
// 2 sqrt(2) 2^-(factor_bits - 1) = 2^-122.5.
//
// A factor between longer parts. A part is m 2^q with m of at most 128 bits, held with its leading bit at bit
// 127. Where the two exponents are within long_difference_shift the difference is exact as a signed 256-bit
// integer (m 2^126 and the other part: below 2^255). Where they are farther apart the part of smaller
// exponent is below 2^(q_large + 2), a quarter of the held factor's unit (the larger part alone has its
// leading bit at q_large + 127, and the unit is factor_bits below a leading bit), and the factor is the larger
// part. Either way the factor is then held as above, each part within two units of the exact one, and the
// bound below is the same.
//
// The product. It is held as two signed 128-bit integers and an exponent, the larger part in
// [2^(product_bits - 1), 2^product_bits). A step forms the exact complex product with a factor in 256 bits
// (each part below 2^(product_bits + factor_bits + 1)) and shifts it back to product_bits, a floor of each part:
// below one unit each where the modulus is at least 2^(product_bits - 1), relative error at most
// sqrt(2) 2^-(product_bits - 1) = 2^-124.5.
//
// The bound. With rho the relative error of the held product (against the modulus of the exact one), a step
// gives (1 + rho)(1 + 2^-122.5)(1 + 2^-124.5) - 1, so after m factors rho <= (1 + 2^-122.18)^m - 1, which is below
// m 2^step_error_log2 for every m below 2^60, and the exact product is within rho / (1 - rho) <= m 2^-121 of the
// held one in units of the held one's modulus (the factor two covers 2^-122.18 against 2^step_error_log2 = 2^-122
// and the division). product() returns the held value with that radius in both coordinates.
class ShortNodeProducts {
  public:
    // A binary64 number as mantissa 2^exponent with 2^52 <= |mantissa| < 2^53, or zero (mantissa zero).
    struct Part {
        int64_t mantissa = 0;
        int exponent = 0;
    };
    // A number of at most two limbs as magnitude 2^exponent, the magnitude's leading bit at bit 127, or zero.
    struct LongPart {
        uint128 magnitude = 0;
        bool negative = false;
        slong exponent = 0;
    };

  private:
    // Every node's parts as long parts; and as binary64 parts where every part of every node is one.
    std::vector<LongPart> long_real_, long_imaginary_;
    std::vector<Part> real_, imaginary_;
    bool usable_ = false, binary64_ = false;

  public:
    // Bits a binary64 part is shifted up by at most when the two exponents of a difference are far apart; bits
    // of a held factor's larger part; of the held product's larger part. A difference (53 + difference_shift
    // bits and one for the sum), factor_bits + 1 and product_bits + 1 fit 127 bits, and a step's product
    // (product_bits + factor_bits + 2 bits) fits 255.
    static constexpr int difference_shift = 73, factor_bits = 125, product_bits = 126;
    // The same for longer parts: 128 + long_difference_shift bits and one for the sum fit 255.
    static constexpr int long_difference_shift = 126;
    // A step's relative error is at most 2^-122.5 + 2^-124.5 + their product, below 2^step_error_log2.
    static constexpr int step_error_log2 = -122;
    // The largest exponent of two a node's nonzero part may have in size (and the least, negated): half the
    // binary64 range, so that a binary64 part is a normal number with room to spare and is read exactly.
    static constexpr int exponent_range = 500;

    // The nodes must be exact balls. usable(): every part of every node is zero or a number of at most two
    // limbs within exponent_range, and there are at least two nodes. binary64(): every part is a binary64
    // number, and the factors take the shorter form.
    explicit ShortNodeProducts(std::span<const ComplexBall> nodes);
    bool usable() const { return usable_; }
    bool binary64() const { return binary64_; }
    // The relative accuracy, in bits, of a product of n nodes (n - 1 factors): -log2((n - 1) 2^(step_error_log2 + 1)),
    // rounded down.
    static slong accuracy_bits(size_t nodes);
    // prod_{j != i} (b_i - b_j) as a ball, and a lower bound of the least |b_i - b_j|. False when node i
    // coincides with another one (a zero factor), as a ball product whose enclosure contains zero is refused.
    bool product(acb_t product, mag_t separation, size_t i) const;
    // prod over j in `subset` of (b_i - c_j), c_j node j of `others` (another usable set of points, or this
    // one): the factors by which a model's residue at an unmoved node changes when the nodes of `subset` move.
    // Its radius is the bound's for subset.size() factors. False for a zero factor or an index outside
    // `others`.
    bool product_over(acb_t product, size_t i, std::span<const size_t> subset, const ShortNodeProducts &others) const;
};
} // namespace certpolysolve::detail
