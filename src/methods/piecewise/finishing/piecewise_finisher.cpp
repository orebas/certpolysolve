// PiecewiseRootState::Impl's source finisher (the class: piecewise_state.hpp): at a batch boundary, owned roots
// refined by Newton steps on the input itself (begin_finishing, advance_finisher, observe_finisher), each finished
// owner marked (mark_finished_owner). Design: docs/piecewise-root-certification.md (the finisher).
#include "methods/piecewise/piecewise_state.hpp"
#include "methods/piecewise/finishing/piecewise_finisher_source.hpp"
#include "core/resource_limits.hpp"
#include <cmath>


namespace certpolysolve::detail {

bool PiecewiseRootState::Impl::owner_finished(size_t id) const {
    return id < accurate_owner.size() && accurate_owner[id];
}
void PiecewiseRootState::Impl::mark_finished_owner(const Root &proof) {
    if (!config.source_finisher || !ownership || proof.count != 1) return;
    StageTimer timer(stats.finisher_ns);
    const auto id = ownership->identifying(proof);
    if (!id) return;
    accurate_owner.resize(std::max(accurate_owner.size(), *id + 1), false);
    accurate_owner[*id] = true;
    unfinished_owners.erase(*id);
}
void PiecewiseRootState::Impl::begin_finishing(double next_work, const cps_options &o) {
    if (!config.source_finisher || !ownership || !(next_work > 0) ||
        o.cancelled() || finishing || !covered || supplied) return;
    StageTimer total_timer(stats.finisher_ns);
    StageTimer admission_timer(stats.finisher_admission_ns);
    const PiecewiseWork work_before = stats.work;
    struct Ledger {
        PiecewiseRootStatistics &stats; const PiecewiseWork &before;
        ~Ledger() { PiecewiseWork spent = stats.work; spent -= before; stats.finishing_work += spent; }
    } ledger{stats, work_before};
    if (failed_finisher_precision == precision && failed_finisher_covered == covered &&
        failed_finisher_improvements == stats.candidate_improvements) return;
    if (finisher_precision != precision) {
        finisher_precision = precision;
        finisher_source.reset();
        pilot_work = marginal_pilot_work = 0;
    }
    Finishing task;
    task.budget = next_work;
    const double walk_before = stats.work.walk;
    for (auto it = unfinished_owners.begin(); it != unfinished_owners.end();) {
        const size_t id = *it;
        stats.work.walk += 1;
        if (!ownership->active(id) || owner_finished(id)) { it = unfinished_owners.erase(it); continue; }
        ++it;
        const size_t candidate = owner_candidate[id];
        if (candidate < carried.size() && owner_of(candidate) == id &&
            ownership->contains(id, carried[candidate].value())) task.owners.emplace_back(id, candidate);
    }
    if (task.owners.empty()) return;
    // Source preparation is charged once to the pilot's batch, not once
    // per root in its extrapolation. Check its structural
    // storage estimate before copying the owned coefficient representation.
    // Bound the source snapshot, block polynomial/derivative, powers and
    // cached second-derivative magnitudes conservatively per coefficient.
    const size_t complex_values = 4;
    const size_t mantissas = 12;
    const auto coefficient_bytes = resources::add(sizeof(Coefficient) + complex_values * sizeof(acb_struct),
        resources::multiply(static_cast<uint64_t>(precision) / 8 + 1, mantissas));
    if (!resources::storage_fits(resources::multiply(length, coefficient_bytes), o.temporary_bytes)) {
        ++stats.finisher_refusals;
        return;
    }
    task.spent = task.admission = piecewise_work_factors.walk * (stats.work.walk - walk_before);
    if (task.spent >= task.budget) { ++stats.finisher_refusals; return; }
    if (marginal_pilot_work > 0 &&
        task.spent + marginal_pilot_work * static_cast<double>(task.owners.size()) > task.budget) {
        ++stats.finisher_refusals;
        observe_finisher(task, "refused-estimate");
        return;
    }
    finishing = std::move(task);
    ++stats.finisher_attempts;
    observe_finisher(*finishing, "admitted");
}
void PiecewiseRootState::Impl::observe_finisher(const Finishing &task, const char *phase, double step, bool certificate) const {
    if (!config.observe_finisher) return;
    PiecewiseFinisherEvent event;
    event.source = config.diagnostic_source_id;
    event.arithmetic_bits = static_cast<uint64_t>(precision);
    event.covered = covered;
    event.owners = task.owners.size();
    event.position = task.position;
    event.budget = task.budget;
    event.spent = task.spent;
    event.pilot = pilot_work;
    event.preparation = task.preparation;
    event.admission = task.admission;
    event.marginal_pilot = marginal_pilot_work;
    event.step = step;
    event.certificate = certificate;
    config.observe_finisher(event, phase);
}
bool PiecewiseRootState::Impl::admit_source_certificate(PiecewiseSourceCertificate certificate, size_t owner) {
    if (!finisher_source || !certificate.belongs_to(*finisher_source))
        throw std::logic_error("piecewise certificate belongs to another source");
    const Root &proof = certificate.proof();
    if (proof.count != 1 || !ownership->identifies(owner, proof) || index.overlaps(proof, result.roots))
        return false;
    if (covered >= static_cast<uint64_t>(result.degree))
        throw std::logic_error("original-source finisher exceeds input degree");
    index.insert(proof, result.roots.size());
    result.roots.push_back(std::move(certificate).release());
    ++covered;
    ++stats.certificates;
    ++stats.finisher_certificates;
    accurate_owner[owner] = true;
    unfinished_owners.erase(owner);
    if (!guards.insert(ownership->proof(owner), result.roots.back()))
        guards.insert(result.roots.back(), result.roots.back());
    publish_reflection();
    return true;
}
void PiecewiseRootState::Impl::advance_finisher(const cps_options &o) {
    if (!finishing || o.cancelled()) return;
    StageTimer step_timer(stats.finisher_ns);
    const PiecewiseWork ledger_before = stats.work;
    struct Ledger {
        PiecewiseRootStatistics &stats; const PiecewiseWork &before;
        ~Ledger() { PiecewiseWork spent = stats.work; spent -= before; stats.finishing_work += spent; }
    } ledger{stats, ledger_before};
    const double work_before = counted_work(stats.work);
    const double preparation_before = stats.work.finisher_preparation;
    auto &task = *finishing;
    bool stop = false, no_progress = false, new_certificate = false;
    const auto step = [&] {
        if (task.position == task.owners.size()) { stop = true; return; }
        const auto [owner, candidate] = task.owners[task.position];
        if (!ownership->active(owner) || owner_of(candidate) != owner || owner_finished(owner)) {
            ++task.position;
            return;
        }
        if (!finisher_source) {
            StageTimer preparation(stats.finisher_preparation_ns);
            finisher_source = std::make_unique<PiecewiseFinisherSource>(
                std::vector<Coefficient>(input.coefficients.begin(),
                    input.coefficients.begin() + static_cast<ptrdiff_t>(length)), stats);
        }
        ComplexBall value, derivative, correction, proposal;
        if (!finisher_source->evaluate(value.value(), derivative.value(), carried[candidate].value(), precision) ||
            acb_contains_zero(derivative.value())) { no_progress = true; return; }
        if (o.cancelled()) return;
        acb_div(correction.value(), value.value(), derivative.value(), precision);
        acb_sub(proposal.value(), carried[candidate].value(), correction.value(), precision);
        acb_get_mid(proposal.value(), proposal.value());
        if (!acb_is_finite(proposal.value()) || !ownership->contains(owner, proposal.value())) {
            no_progress = true;
            return;
        }
        const bool moved = !acb_equal(proposal.value(), carried[candidate].value());
        carried[candidate] = proposal;
        auto certificate = finisher_source->certify(proposal.value(),o,precision);
        if (o.cancelled()) return;
        if (certificate && admit_source_certificate(std::move(*certificate), owner)) {
            ++task.position;
            task.last_correction = std::numeric_limits<double>::infinity();
            new_certificate = true;
            return;
        }
        Magnitude magnitude;
        acb_get_mag(magnitude.value(), correction.value());
        const double level = std::ceil(mag_get_d_log2_approx(magnitude.value()));
        no_progress = !moved || !(level < task.last_correction);
        task.last_correction = level;
    };
    step();
    const double spent = counted_work(stats.work) - work_before;
    task.spent += spent;
    task.preparation += piecewise_work_factors.finisher_preparation *
                        (stats.work.finisher_preparation - preparation_before);
    if (new_certificate && !(pilot_work > 0)) {
        pilot_work = task.spent;
        const double setup = task.admission + task.preparation;
        // The total pilot remains observable. Only its per-root work is
        // extrapolated; all setup still consumes the actual batch budget.
        marginal_pilot_work = pilot_work > setup ? pilot_work - setup : pilot_work;
        // Keep the pilot's valid certificate even if the extrapolated
        // remaining work loses against completing the suspended model batch.
        if (marginal_pilot_work * static_cast<double>(task.owners.size() - task.position) + task.spent >
            task.budget) { stop = true; ++stats.finisher_refusals; }
        observe_finisher(task, "pilot", spent, new_certificate);
    }
    if (no_progress) {
        failed_finisher_precision = precision;
        failed_finisher_covered = covered;
        failed_finisher_improvements = stats.candidate_improvements;
    }
    if (no_progress || stop || task.position == task.owners.size() || task.spent >= task.budget) {
        observe_finisher(task, no_progress ? "no-progress" : stop ? "stopped" :
            task.position == task.owners.size() ? "cohort-finished" : "budget-spent", spent, new_certificate);
        finishing.reset();
    }
}

} // namespace certpolysolve::detail
