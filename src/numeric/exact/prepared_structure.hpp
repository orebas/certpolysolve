#pragma once
#include "numeric/exact/input_structure.hpp"
#include <vector>

namespace certpolysolve::detail {
enum class PreparedRepresentation { Compressed, Integer, Gaussian };

// A value-only description of the *actual* retained numerical factor. Shape
// may reveal further possible reductions; power records only an applied one.
// No coefficient/proof objects or references to a mutable solver are retained.
struct PreparedFactorStructure {
    InputStructure shape;
    uint64_t multiplicity = 1, power = 1;
    bool squarefree_proved = false;
    PreparedRepresentation representation = PreparedRepresentation::Compressed;
};

struct PreparedStructure {
    int64_t input_degree = -1;
    uint64_t zero_order = 0;
    std::vector<PreparedFactorStructure> factors;

    // p(x) = c x^s product_j q_j(x^k_j)^m_j. Each shape describes q_j
    // after its exact zero factor, if any, is accounted for separately.
    // Budget-limited raw factors remain explicitly unproved square-free.
    void append(PreparedFactorStructure factor) {
        const auto weight = resources::multiply(factor.power, factor.multiplicity);
        zero_order = resources::add(zero_order, resources::multiply(factor.shape.start, weight));
        if (factor.shape.degree()) factors.push_back(std::move(factor));
    }
    uint64_t accounted_degree() const {
        uint64_t degree = zero_order;
        for (const auto &factor : factors)
            degree = resources::add(degree, resources::multiply(factor.shape.degree(),
                resources::multiply(factor.power, factor.multiplicity)));
        return degree;
    }
    uint64_t numerical_degree() const {
        uint64_t degree = 0;
        for (const auto &factor : factors) degree = resources::add(degree, factor.shape.degree());
        return degree;
    }
};

InputStructure input_structure(const IntegerPolynomial &);
InputStructure input_structure(std::span<const GaussianRational>);
} // namespace certpolysolve::detail
