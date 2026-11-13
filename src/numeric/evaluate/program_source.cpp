// The polynomial a program computes, as the MonicSource the structured route iterates on and certifies with
// (ProgramSource, declared in program_source.hpp): its value and derivative over a ball, in complex balls at the
// precision asked.
//
// ProgramSource::evaluate runs the program once over DualBallArithmetic through run_prepared (the last uses computed
// once, at construction, and the values' storage kept between calls), then divides the value and the derivative by the
// exact leading coefficient the preparation proved (skipped when it is 1; its ball kept at the last precision). With a
// zero order v > 0 the source is q = p / z^v, the program's polynomial without its root of order v at zero: q = p / z^v
// and q' = (p' - v p / z) / z^v, which are not finite on a ball that holds zero. A value or a derivative that is not
// finite makes evaluate return false. Nothing else is known about the polynomial: taylor() and the coefficient views
// refuse (MonicSource's defaults), and real() is false whatever the constants.
//
// State: the program, its prepared form and the values' storage, and the leading coefficient's ball at the last
// precision.
// Entry points: the constructor and evaluate, from solver_program.cpp (program_advance: the source of the preliminary
// phase and the refinement, the Newton-disk pass's source and, for a root at zero, its undivided p) and
// program_recurrence_seeds.cpp (the circle counts of a sub-program no structural rule seeds).
// Not here: the batched machine evaluator of the proposals (program_machine.cpp), the preparation that proves the
// degree and leading coefficient (program_preparation.cpp).
// Design notes: docs/slp-input-spec-20261015.md (section 6, the solve), docs/slp-second-version-note-20261017.md (the
// Newton-disk pass that certifies on this source).
#include "numeric/evaluate/program_arithmetic.hpp"
#include "numeric/evaluate/program_interpreter.hpp"
#include "numeric/evaluate/program_source.hpp"
#include <stdexcept>

namespace certpolysolve::detail {

ProgramSource::ProgramSource(Program program, slong degree, GaussianRational leading, uint64_t zero_order)
    : program_(std::move(program)), prepared_(program_), degree_(degree), leading_(std::move(leading)),
      leading_is_one_(fmpq_is_one(leading_.real.value()) && fmpq_is_zero(leading_.imag.value())),
      zero_order_(zero_order) {
    if (degree_ < 1 || leading_.zero())
        throw std::invalid_argument("a program source needs a proved degree of at least 1 and a nonzero leading coefficient");
}

bool ProgramSource::evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) {
    DualBallArithmetic arithmetic(z, precision);
    auto result = run_prepared(prepared_, arithmetic, values_);
    if (leading_is_one_) {
        acb_swap(value, result.value.value());
        acb_swap(derivative, result.derivative.value());
    } else {
        if (leading_precision_ != precision) {
            arb_set_fmpq(acb_realref(leading_ball_.value()), leading_.real.value(), precision);
            arb_set_fmpq(acb_imagref(leading_ball_.value()), leading_.imag.value(), precision);
            leading_precision_ = precision;
        }
        acb_div(value, result.value.value(), leading_ball_.value(), precision);
        acb_div(derivative, result.derivative.value(), leading_ball_.value(), precision);
    }
    if (zero_order_) {
        // q = p / z^v and q' = (p' - v p / z) / z^v: p without its root of order v at zero.
        ComplexBall power, quotient;
        acb_pow_ui(power.value(), z, zero_order_, precision);
        acb_div(quotient.value(), value, z, precision);
        acb_mul_ui(quotient.value(), quotient.value(), zero_order_, precision);
        acb_sub(derivative, derivative, quotient.value(), precision);
        acb_div(value, value, power.value(), precision);
        acb_div(derivative, derivative, power.value(), precision);
    }
    return acb_is_finite(value) && acb_is_finite(derivative);
}

} // namespace certpolysolve::detail
