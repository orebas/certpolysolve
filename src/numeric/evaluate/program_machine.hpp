#pragma once
// A program (an SLP, program.hpp) evaluated in machine arithmetic at many points at once, for the proposals of
// the program route's binary64 stage and its root counts (program_seeds): value and derivative by forward-mode
// differentiation, and a noise estimate per point. Nothing here certifies: the balls of ProgramSource do.
//
// Why batched: the generic interpreter (run_prepared) evaluates one point at a time through a switch per node
// with optional slots, several times slower than a hand-written loop; here every node is one loop over a chunk
// of points, so the dispatch is paid once a node and the loops vectorize. Why its own copy of the program: the
// output's dependency cone only (a dead node's overflow would otherwise poison a point for nothing), the
// constants rounded once per arithmetic, no borrowed Program.
// History: docs/slp-seed-placement-note-20261016.md, section 9.
//
// The noise estimate, the floor a point's evaluation cannot be read below (binary64 stage: a value within a
// few times its noise is at the roundoff floor, and the point is handed to the balls). It is the first-order
// propagated rounding of the LAST sum on each path: a sum or difference starts afresh at u (|a| + |b|), its
// operands' noise discarded; a product propagates (|a| n_b + |b| n_a, plus its own rounding); a power is its
// products. This is a lower bound of the error, by design. A bound that propagates through sums grows like the
// balls' radii at a recurrence's roots (about 1.27 n bits for Chebyshev T_n at |z| = 1: the radii follow the
// recurrence with |2z| + 1 while the values stay O(1)) and would call every point of every recurrence unreadable
// in doubles. What the lower bound catches is cancellation at the output's last sums (Laguerre's
// ((2k+1-x) L_k - k L_(k-1)) / (k+1) at its large roots: u (|A| + |B|) against a value far below it); what it
// misses costs binary64 sweeps only, bounded by the stage's window rule and the handover read in balls
// (program_advance, solver_program.cpp).
//
// Same bits on every CPU (a project rule): the arithmetic loops are defined in program_machine.cpp only, under the library's
// -ffp-contract=off, over machine::Complex<R> (a fixed order of operations), with no std::complex and no fma;
// every loop is elementwise (element j depends on element j only), so the compiler's vectorization keeps each
// element's order. evaluate_one goes through the same loops with one point: the reference the batched loops
// equal bit for bit (tests/test_program.cpp, machine_bits).
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "numeric/arithmetic/block_double.hpp"
#include "input/program/program.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace certpolysolve::detail {

template <class R> class ProgramMachine {
  public:
    using C = machine::Complex<R>;
    struct Evaluation {
        C value, derivative;
        R noise, noise_derivative; // the floor of |value| and |derivative|; 0 for a constant program
        bool finite = false;       // value and derivative finite (a program has no division: inf and NaN
                                   // reach the output)
    };
    // Copies the output's dependency cone of `program` and rounds its constants to R (to the nearest, from a
    // 53-bit enclosure). valid() is false when the program has no output or a constant is outside R's range
    // (a nonzero constant rounding to zero included); nothing else may be called then.
    //
    // With `zero_order` v > 0 the evaluator stands for q = p / x^v, the program without its root of order v at
    // zero (ProgramSource's zero order), scaled by z^v: the value is p(z) and the derivative p'(z) - v p(z) / z,
    // which are z^v q(z) and z^v q'(z). Every reading of the evaluations is free of that scale: the Newton ratio
    // value / derivative is q / q', the Cauchy term z derivative / value is z q' / q, and each floor test weighs a
    // value against its own noise. z^v itself, which leaves the range where p does not, is never formed. The
    // derivative's noise is p''s, plus v times p's noise over |z|, plus the subtraction's rounding. A point at 0 is
    // not finite.
    explicit ProgramMachine(const Program &program, uint64_t zero_order = 0);
    // The same cone in another arithmetic (the range tier, machine::BlockDouble, from a valid binary64 evaluator:
    // its constants are doubles, which the other arithmetic holds exactly).
    template <class Other> explicit ProgramMachine(const ProgramMachine<Other> &other);
    bool valid() const { return valid_; }
    // The cone's operations: the cost unit of an evaluation (statistics.evaluation_terms).
    size_t terms() const { return nodes_.size(); }
    // The arrays live at once (one set per slot), and the storage of a chunk of `points`.
    size_t slots() const { return slots_; }
    uint64_t storage_bytes(size_t points) const;
    // The program at every point of z, into out (out.size() >= z.size()), in chunks of at most `chunk` points
    // (at least 1); the chunk's arrays are kept between calls, sized to the largest chunk seen.
    // Keep this thin dispatcher visible to callers across module boundaries;
    // the numerical kernel and its fixed arithmetic order remain out of line.
    void evaluate(std::span<const C> z, std::span<Evaluation> out, size_t chunk) {
        chunk = std::max<size_t>(chunk, 1);
        for (size_t start = 0; start < z.size(); start += chunk) {
            const size_t count = std::min(chunk, z.size() - start);
            run_chunk(z.subspan(start, count), out.subspan(start, count));
        }
    }
    // One point through the same loops (a chunk of one).
    Evaluation evaluate_one(C z);

  private:
    struct Arrays {
        enum Field { ValueRe, ValueIm, DerivativeRe, DerivativeIm, Noise, NoiseDerivative, Fields };
        std::vector<R> values;
        size_t stride = 0;
        void resize(size_t points);
        R *field(Field f) { return values.data() + static_cast<size_t>(f) * stride; }
        const R *field(Field f) const { return values.data() + static_cast<size_t>(f) * stride; }
    };
    void run_chunk(std::span<const C> z, std::span<Evaluation> out);
    std::vector<ProgramNode> nodes_; // the cone, renumbered; a Constant's `value` indexes constants_
    std::vector<C> constants_;
    uint64_t output_ = 0;
    std::vector<uint32_t> slot_; // each node's slot
    size_t slots_ = 0;
    bool valid_ = false;
    std::vector<Arrays> arrays_;
    size_t capacity_ = 0;
    uint64_t zero_order_ = 0; // v of q = p / x^v (0: the program itself)
    template <class> friend class ProgramMachine;
};

extern template class ProgramMachine<double>;
extern template class ProgramMachine<machine::BlockDouble>;

} // namespace certpolysolve::detail
