// Batched program values, derivatives and noise estimates, with identical
// per-point operations in every batch size. Candidate iteration is implemented
// separately in numeric/candidates/program_candidates.cpp. No result here certifies.
#include "numeric/evaluate/program_machine.hpp"
#include "numeric/evaluate/program_slots.hpp"
#include "numeric/arithmetic/machine_conversion.hpp"
#include "core/raii.hpp"
#include "core/resource_limits.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace certpolysolve::detail {
namespace {
inline double absolute(double x) { return std::abs(x); }
inline machine::BlockDouble absolute(machine::BlockDouble x) { return machine::abs(x); }
// |re| + |im|: a magnitude without a square root, within a factor sqrt 2 of the modulus, for the noise.
template <class R> inline R magnitude(R re, R im) { return absolute(re) + absolute(im); }

// An odd number of 64-byte lines separates the six channels of a live slot.
// Power-of-two byte strides make these simultaneously accessed arrays compete
// for the same cache sets. This is storage padding only; arithmetic sees count
// points. Tiny batches retain their exact size. No CPU-dependent arithmetic.
template <class R> uint64_t channel_stride(size_t points) {
    constexpr size_t line_bytes = 64;
    static_assert(sizeof(R) <= line_bytes && line_bytes % sizeof(R) == 0);
    constexpr size_t per_line = line_bytes / sizeof(R);
    if (points < per_line)
        return points;
    const auto lines = static_cast<uint64_t>(points / per_line) + (points % per_line != 0);
    return resources::multiply(lines | uint64_t{1}, per_line);
}
} // namespace

template <class R> ProgramMachine<R>::ProgramMachine(const Program &program, uint64_t zero_order) : zero_order_(zero_order) {
    if (!program.output || *program.output >= program.nodes.size())
        return;
    // The output's cone: reachable backwards. Operands are earlier nodes, so one pass downward finds it.
    std::vector<uint8_t> reachable(program.nodes.size(), 0);
    reachable[*program.output] = 1;
    for (size_t i = program.nodes.size(); i-- > 0;) {
        if (!reachable[i])
            continue;
        const auto &node = program.nodes[i];
        if (node.op == ProgramOp::Variable || node.op == ProgramOp::Constant)
            continue;
        reachable[node.a] = 1;
        if (node.op != ProgramOp::Negate && node.op != ProgramOp::Power)
            reachable[node.b] = 1;
    }
    std::vector<uint64_t> renumbered(program.nodes.size(), 0);
    for (size_t i = 0; i < program.nodes.size(); ++i) {
        if (!reachable[i])
            continue;
        renumbered[i] = nodes_.size();
        ProgramNode node = program.nodes[i];
        if (node.op == ProgramOp::Constant) {
            ComplexBall enclosure;
            program.constants[node.value].enclose(enclosure.value(), 53);
            C constant;
            if (!from_ball(constant, enclosure.value()))
                return; // outside R's range: valid_ stays false
            node.value = constants_.size();
            constants_.push_back(constant);
        } else if (node.op != ProgramOp::Variable) {
            node.a = renumbered[node.a];
            if (node.op != ProgramOp::Negate && node.op != ProgramOp::Power)
                node.b = renumbered[node.b];
        }
        nodes_.push_back(node);
    }
    output_ = renumbered[*program.output];
    // The slot plan: each node takes a slot from the free list; an operand's slot is released only after its
    // reader's slot is assigned, so no node writes an operand's arrays.
    Program cone;
    cone.nodes = nodes_;
    cone.output = output_;
    ProgramSlots schedule(cone);
    // Keep the machine schedule's original compact representation. The shared
    // planner also serves the ball interpreter, whose indices are size_t.
    if (schedule.count - 1 > std::numeric_limits<uint32_t>::max())
        return;
    slot_.assign(schedule.slot.begin(), schedule.slot.end());
    slots_ = schedule.count;
    valid_ = true;
}

template <class R>
template <class Other>
ProgramMachine<R>::ProgramMachine(const ProgramMachine<Other> &other)
    : nodes_(other.nodes_), output_(other.output_), slot_(other.slot_), slots_(other.slots_), valid_(other.valid_),
      zero_order_(other.zero_order_) {
    constants_.reserve(other.constants_.size());
    for (const auto &c : other.constants_)
        constants_.push_back(C{R(c.re), R(c.im)});
}

template <class R> uint64_t ProgramMachine<R>::storage_bytes(size_t points) const {
    return resources::multiply(resources::multiply(static_cast<uint64_t>(slots_), channel_stride<R>(points)),
                               static_cast<uint64_t>(Arrays::Fields * sizeof(R)));
}

template <class R> void ProgramMachine<R>::Arrays::resize(size_t points) {
    const auto next = channel_stride<R>(points);
    const auto bytes = resources::multiply(next, Fields * sizeof(R));
    if (!resources::storage_fits(bytes, resources::unlimited))
        throw std::length_error("program machine buffers exceed addressable storage");
    values.resize(static_cast<size_t>(bytes / sizeof(R)));
    stride = static_cast<size_t>(next);
}

template <class R> typename ProgramMachine<R>::Evaluation ProgramMachine<R>::evaluate_one(C z) {
    Evaluation result;
    run_chunk(std::span<const C>(&z, 1), std::span<Evaluation>(&result, 1));
    return result;
}

// Every loop below is elementwise: element j of a node's arrays depends on element j of its operands' arrays
// only, in a fixed order of operations (Complex<R>'s product order for products), so a chunk of one point and
// a chunk of many return the same bits for the same point.
template <class R> void ProgramMachine<R>::run_chunk(std::span<const C> z, std::span<Evaluation> out) {
    const size_t count = z.size();
    if (count > capacity_) {
        arrays_.resize(slots_);
        for (auto &a : arrays_)
            a.resize(count);
        capacity_ = count;
    }
    const R u = R(0x1p-53), one = R(1), zero = R{};
    for (size_t i = 0; i < nodes_.size(); ++i) {
        const auto &node = nodes_[i];
        Arrays &r = arrays_[slot_[i]];
        R *const vre = r.field(Arrays::ValueRe), *const vim = r.field(Arrays::ValueIm);
        R *const dre = r.field(Arrays::DerivativeRe), *const dim = r.field(Arrays::DerivativeIm);
        R *const nv = r.field(Arrays::Noise), *const nd = r.field(Arrays::NoiseDerivative);
        switch (node.op) {
        case ProgramOp::Variable:
            for (size_t j = 0; j < count; ++j) {
                vre[j] = z[j].re;
                vim[j] = z[j].im;
                dre[j] = one;
                dim[j] = zero;
                nv[j] = zero;
                nd[j] = zero;
            }
            break;
        case ProgramOp::Constant: {
            const C c = constants_[node.value];
            for (size_t j = 0; j < count; ++j) {
                vre[j] = c.re;
                vim[j] = c.im;
                dre[j] = zero;
                dim[j] = zero;
                nv[j] = zero;
                nd[j] = zero;
            }
            break;
        }
        case ProgramOp::Add:
        case ProgramOp::Subtract: {
            const Arrays &a = arrays_[slot_[node.a]], &b = arrays_[slot_[node.b]];
            const R *const are = a.field(Arrays::ValueRe), *const aim = a.field(Arrays::ValueIm);
            const R *const adre = a.field(Arrays::DerivativeRe), *const adim = a.field(Arrays::DerivativeIm);
            const R *const bre = b.field(Arrays::ValueRe), *const bim = b.field(Arrays::ValueIm);
            const R *const bdre = b.field(Arrays::DerivativeRe), *const bdim = b.field(Arrays::DerivativeIm);
            if (node.op == ProgramOp::Add) {
                for (size_t j = 0; j < count; ++j) {
                    vre[j] = are[j] + bre[j];
                    vim[j] = aim[j] + bim[j];
                    dre[j] = adre[j] + bdre[j];
                    dim[j] = adim[j] + bdim[j];
                    nv[j] = u * (magnitude(are[j], aim[j]) + magnitude(bre[j], bim[j]));
                    nd[j] = u * (magnitude(adre[j], adim[j]) + magnitude(bdre[j], bdim[j]));
                }
            } else {
                for (size_t j = 0; j < count; ++j) {
                    vre[j] = are[j] - bre[j];
                    vim[j] = aim[j] - bim[j];
                    dre[j] = adre[j] - bdre[j];
                    dim[j] = adim[j] - bdim[j];
                    nv[j] = u * (magnitude(are[j], aim[j]) + magnitude(bre[j], bim[j]));
                    nd[j] = u * (magnitude(adre[j], adim[j]) + magnitude(bdre[j], bdim[j]));
                }
            }
            break;
        }
        case ProgramOp::Multiply: {
            const Arrays &a = arrays_[slot_[node.a]], &b = arrays_[slot_[node.b]];
            const R *const are = a.field(Arrays::ValueRe), *const aim = a.field(Arrays::ValueIm);
            const R *const adre = a.field(Arrays::DerivativeRe), *const adim = a.field(Arrays::DerivativeIm);
            const R *const an = a.field(Arrays::Noise), *const and_ = a.field(Arrays::NoiseDerivative);
            const R *const bre = b.field(Arrays::ValueRe), *const bim = b.field(Arrays::ValueIm);
            const R *const bdre = b.field(Arrays::DerivativeRe), *const bdim = b.field(Arrays::DerivativeIm);
            const R *const bn = b.field(Arrays::Noise), *const bnd = b.field(Arrays::NoiseDerivative);
            for (size_t j = 0; j < count; ++j) {
                // a b, then a' b + a b' (the product rule), each product in Complex<R>'s order.
                const R pre = are[j] * bre[j] - aim[j] * bim[j], pim = are[j] * bim[j] + aim[j] * bre[j];
                const R qre = adre[j] * bre[j] - adim[j] * bim[j], qim = adre[j] * bim[j] + adim[j] * bre[j];
                const R rre = are[j] * bdre[j] - aim[j] * bdim[j], rim = are[j] * bdim[j] + aim[j] * bdre[j];
                vre[j] = pre;
                vim[j] = pim;
                dre[j] = qre + rre;
                dim[j] = qim + rim;
                const R ma = magnitude(are[j], aim[j]), mb = magnitude(bre[j], bim[j]);
                const R mda = magnitude(adre[j], adim[j]), mdb = magnitude(bdre[j], bdim[j]);
                nv[j] = ma * bn[j] + mb * an[j] + u * magnitude(pre, pim);
                nd[j] = mda * bn[j] + mb * and_[j] + ma * bnd[j] + mdb * an[j] +
                        u * (magnitude(qre, qim) + magnitude(rre, rim));
            }
            break;
        }
        case ProgramOp::Negate: {
            const Arrays &a = arrays_[slot_[node.a]];
            const R *const are = a.field(Arrays::ValueRe), *const aim = a.field(Arrays::ValueIm);
            const R *const adre = a.field(Arrays::DerivativeRe), *const adim = a.field(Arrays::DerivativeIm);
            const R *const an = a.field(Arrays::Noise), *const and_ = a.field(Arrays::NoiseDerivative);
            for (size_t j = 0; j < count; ++j) {
                vre[j] = -are[j];
                vim[j] = -aim[j];
                dre[j] = -adre[j];
                dim[j] = -adim[j];
                nv[j] = an[j];
                nd[j] = and_[j];
            }
            break;
        }
        case ProgramOp::Power: {
            // a^e = lower a with lower = a^(e-1) by binary powering (no squaring past the last exponent bit);
            // the derivative e lower a'. The noise follows each product taken.
            const Arrays &a = arrays_[slot_[node.a]];
            const R *const are = a.field(Arrays::ValueRe), *const aim = a.field(Arrays::ValueIm);
            const R *const adre = a.field(Arrays::DerivativeRe), *const adim = a.field(Arrays::DerivativeIm);
            const R *const an = a.field(Arrays::Noise), *const and_ = a.field(Arrays::NoiseDerivative);
            const uint64_t exponent = node.value;
            const R e = R(static_cast<double>(exponent));
            for (size_t j = 0; j < count; ++j) {
                R lre = one, lim = zero, ln = zero;                 // lower and its noise
                R bre = are[j], bim = aim[j], bn = an[j];           // the running power of a and its noise
                for (uint64_t bits = exponent - 1; bits;) {
                    if (bits & 1) {
                        const R tre = lre * bre - lim * bim, tim = lre * bim + lim * bre;
                        ln = magnitude(lre, lim) * bn + magnitude(bre, bim) * ln + u * magnitude(tre, tim);
                        lre = tre;
                        lim = tim;
                    }
                    bits >>= 1;
                    if (bits) {
                        const R sre = bre * bre - bim * bim, sim = bre * bim + bim * bre;
                        bn = (magnitude(bre, bim) + magnitude(bre, bim)) * bn + u * magnitude(sre, sim);
                        bre = sre;
                        bim = sim;
                    }
                }
                const R pre = lre * are[j] - lim * aim[j], pim = lre * aim[j] + lim * are[j];
                const R qre = lre * adre[j] - lim * adim[j], qim = lre * adim[j] + lim * adre[j];
                vre[j] = pre;
                vim[j] = pim;
                dre[j] = e * qre;
                dim[j] = e * qim;
                const R ml = magnitude(lre, lim);
                nv[j] = ml * an[j] + magnitude(are[j], aim[j]) * ln + u * magnitude(pre, pim);
                nd[j] = e * (ml * and_[j] + magnitude(adre[j], adim[j]) * ln) + u * magnitude(dre[j], dim[j]);
            }
            break;
        }
        }
    }
    const Arrays &o = arrays_[slot_[output_]];
    const R order = R(static_cast<double>(zero_order_));
    for (size_t j = 0; j < count; ++j) {
        Evaluation &result = out[j];
        result.value = {o.field(Arrays::ValueRe)[j], o.field(Arrays::ValueIm)[j]};
        result.derivative = {o.field(Arrays::DerivativeRe)[j], o.field(Arrays::DerivativeIm)[j]};
        result.noise = o.field(Arrays::Noise)[j];
        result.noise_derivative = o.field(Arrays::NoiseDerivative)[j];
        if (zero_order_) {
            // z^v q and z^v q' for q = p / z^v: p, and p' - v p / z (the constructor's comment).
            const C quotient = result.value / z[j];
            const C term{order * quotient.re, order * quotient.im};
            result.noise_derivative = result.noise_derivative + order * result.noise / magnitude(z[j].re, z[j].im) +
                                      u * (magnitude(result.derivative.re, result.derivative.im) +
                                           magnitude(term.re, term.im));
            result.derivative = result.derivative - term;
        }
        result.finite = machine::finite(result.value) && machine::finite(result.derivative);
    }
}

template class ProgramMachine<double>;
template class ProgramMachine<machine::BlockDouble>;
template ProgramMachine<machine::BlockDouble>::ProgramMachine(const ProgramMachine<double> &);

} // namespace certpolysolve::detail
