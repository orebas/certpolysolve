#pragma once
// Program input, milestone 3 (docs/slp-input-spec-20261015.md, section 6): the monic polynomial a program computes,
// as a MonicSource for the structured route. A value and its derivative come from running the program over dual
// complex balls at the precision asked (DualBallArithmetic) and dividing by the exact leading coefficient the
// preparation proved (section 4). Nothing else is known about the polynomial: taylor() and the coefficient views
// refuse (MonicSource's defaults), so the route's coefficient-only aids switch off (the certificate attempt, the
// preliminary phase, the local and geometric restarts; the unresolved components restart on circles). A program is
// treated as non-real whatever its constants (section 2).
// History: docs/slp-input-spec-20261015.md, section 13.
#include "numeric/evaluate/program_arithmetic.hpp"
#include "numeric/evaluate/program_interpreter.hpp"
#include "input/program/program.hpp"
#include "numeric/evaluate/monic_source.hpp"
#include <optional>
#include <vector>

namespace certpolysolve::detail {

class ProgramSource final : public MonicSource {
    Program program_;
    PreparedProgram prepared_; // the program's last uses, computed once (points at program_)
    slong degree_;
    GaussianRational leading_;
    bool leading_is_one_;
    // The values' storage, kept between evaluations (run_prepared).
    std::vector<std::optional<DualBallArithmetic::Value>> values_;
    // The leading coefficient enclosed at the precision of the last evaluation.
    ComplexBall leading_ball_;
    slong leading_precision_ = 0;
    // The order v of the program's root at zero that the source divides out (program_zero_order; 0: none).
    uint64_t zero_order_;

  public:
    // `degree` and `leading` are the preparation's (prepare_program): the degree is proved and the leading
    // coefficient is nonzero. With `zero_order` v > 0 the source is q = p / x^v, of degree `degree` (p's degree less
    // v; q's leading coefficient is p's): every stage that reads it solves q, and the route publishes the point 0
    // beside q's cover. q(z) and q'(z) are not finite on a ball that holds zero.
    ProgramSource(Program program, slong degree, GaussianRational leading, uint64_t zero_order = 0);
    ProgramSource(const ProgramSource &) = delete;
    ProgramSource &operator=(const ProgramSource &) = delete;
    slong degree() const override { return degree_; }
    bool real() const override { return false; }
    SourceMetadata metadata() const override {
        return {degree_, false, true, SourceRepresentation::Program, InputSemantics::Exact, false};
    }
    // The cost model's unit for a program is its operation count (one dual-ball operation a node).
    uint64_t terms() const override { return static_cast<uint64_t>(program_.nodes.size()); }
    bool evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) override;
    const Program &program() const { return program_; }
    // The order v divided out (0: none): what a machine evaluator of program() divides out to stand for this source.
    uint64_t zero_order() const { return zero_order_; }
};

} // namespace certpolysolve::detail
