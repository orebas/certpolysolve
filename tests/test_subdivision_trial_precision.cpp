#include "support/cost-progress/prepared_isolation.hpp"
#include "reference_math.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <filesystem>
#include <iostream>

namespace d = certpolysolve::detail;
namespace ir = isolation_research;
using reference::Integer;
using reference::Rational;
using reference::require;
struct Point { Rational real, imag; };

namespace {
constexpr uint64_t initial_bits = CPS_DEFAULT_BITS + 32; // maintained coefficient-route start
constexpr uint64_t test_turn_limit = 20000; // a bounded test, not an algorithmic threshold
constexpr double resume_allowance = std::numeric_limits<double>::max() / 4;

Rational endpoint(const d::Endpoint &end) {
    const auto text = end.text();
    return Rational(Integer(text.mantissa)) * reference::power2(std::stoi(text.exponent));
}
std::vector<Point> known_roots(const std::string &name) {
    const Rational a = reference::power2(-128);
    if (name == "symmetric-n2-l128") return {{1-a, 0}, {1+a, 0}};
    if (name == "translated-n2-h180-l128") {
        const Rational center = reference::power2(180);
        return {{center, 0}, {center+a, 0}};
    }
    if (name == "translated-symmetric-n2-h180-l128") {
        const Rational center = reference::power2(180);
        return {{center-a, 0}, {center+a, 0}};
    }
    if (name == "gaussian-n2-l128") return {{-a, 1}, {a, 1}};
    throw std::invalid_argument("unknown exact test fixture");
}
void check_input(const std::filesystem::path &path, const std::vector<Point> &roots) {
    const auto rows = scaling_io::read(path.string());
    require(rows.size() == roots.size()+1, "distinct planted roots must match the degree");
    for (size_t i = 0; i < roots.size(); ++i) {
        for (size_t j = 0; j < i; ++j)
            require(roots[i].real != roots[j].real || roots[i].imag != roots[j].imag, "roots must be distinct");
        Point value{0, 0};
        for (auto coefficient = rows.rbegin(); coefficient != rows.rend(); ++coefficient) {
            const Rational real = value.real*roots[i].real-value.imag*roots[i].imag + Rational(Integer((*coefficient)[0]));
            value.imag = value.real*roots[i].imag+value.imag*roots[i].real + Rational(Integer((*coefficient)[1]));
            value.real = real;
        }
        require(value.real == 0 && value.imag == 0, "independent exact substitution must vanish");
    }
}
void check_cover(const std::vector<d::Root> &cover, const std::vector<Point> &roots) {
    std::vector<unsigned> owned(roots.size(), 0);
    std::vector<std::array<Rational, 4>> rectangles;
    for (const auto &region : cover) {
        std::array<Rational, 4> box;
        for (size_t j=0; j<4; ++j) box[j]=endpoint(region.endpoints[j]);
        require(box[0]<=box[1] && box[2]<=box[3], "exported intervals must be ordered");
        uint64_t count=0;
        for (size_t i=0; i<roots.size(); ++i)
            if (box[0]<=roots[i].real && roots[i].real<=box[1] &&
                box[2]<=roots[i].imag && roots[i].imag<=box[3]) { ++count; ++owned[i]; }
        require(count == region.count && count > 0, "independent planted-root rectangle count differs");
        require(region.multiplicity == (count == 1 ? 1u : 0u), "cluster multiplicity must remain unknown");
        for (const auto &other : rectangles)
            require(box[1]<other[0] || other[1]<box[0] || box[3]<other[2] || other[3]<box[2],
                    "closed exported rectangles must be disjoint");
        rectangles.push_back(std::move(box));
    }
    require(std::all_of(owned.begin(), owned.end(), [](unsigned n) { return n==1; }), "every root must be covered once");
}

struct Observed {
    bool adopted=false, discarded=false, branched=false;
};
cps_stop_reason tick(ir::Run &run, Observed &events, double allowance=resume_allowance) {
    const auto before=run.state();
    const auto previous_enclosures=run.enclosure_calls;
    const auto previous_sum=run.enclosure_bits_sum;
    const auto successes=run.statistics.newton_successes;
    const auto reason=run.tick(allowance);
    const auto after=run.state();
    if (before.size()==1 && !before[0].blocked && run.enclosure_calls>previous_enclosures) {
        const uint64_t work=before[0].trial_bits ? before[0].trial_bits : before[0].base_bits;
        require(run.enclosure_bits_sum-previous_sum==work, "enclosure must use active trial precision");
    }
    if (before.size()==1 && after.size()==1 && before[0].trial_bits>before[0].base_bits) {
        if (run.statistics.newton_successes>successes) {
            require(after[0].base_bits==before[0].trial_bits && !after[0].trial_bits,
                    "successful trial must adopt its actual precision");
            events.adopted=true;
        } else if (after[0].splitting && !after[0].trial_bits) {
            require(after[0].base_bits==before[0].base_bits, "definitive failure must discard trial precision");
            events.discarded=true;
        }
    }
    if (before.size()==1 && before[0].splitting && after.size()>1) {
        for (const auto &child : after)
            require(child.base_bits==before[0].base_bits && !child.trial_bits && child.newton_bits==2,
                    "regrouped children must inherit base precision and reset exponent");
        events.branched=true;
    }
    return reason;
}
void finish(ir::Run &run, Observed &events, const std::vector<Point> &roots) {
    while (!run.goal_met && run.turns<test_turn_limit)
        require(tick(run, events)==CPS_STOP_NONE, "isolation unexpectedly stopped before its goal");
    require(run.goal_met && run.engine().pending()==0, "explicit service must retire at first singleton certificates");
    check_cover(run.cover, roots);
    for (const auto &state : run.state())
        require(state.finished && state.count==1 && !state.trial_bits, "finished roots retain no private trial");
}
void full_isolation(const std::filesystem::path &directory, const std::string &name) {
    const auto path=directory/(name+".ps");
    const auto roots=known_roots(name);
    check_input(path, roots);
    ir::Run run(ir::IntegerFactor(path.string()), initial_bits, ir::isolation_options(4096), d::SubdivisionPrecision::TrialLocal);
    check_cover(run.cover, roots);
    Observed events;
    finish(run, events, roots);
    std::cout << name << " goal=1 turns=" << run.turns << " peak=" << run.peak_enclosure_bits
              << " adopted=" << events.adopted << " discarded=" << events.discarded << '\n';
}
void cap_and_suspension(const std::filesystem::path &directory) {
    const auto path=directory/"translated-symmetric-n2-h180-l128.ps";
    const auto roots=known_roots("translated-symmetric-n2-h180-l128");
    check_input(path, roots);
    ir::Run run(ir::IntegerFactor(path.string()), initial_bits, ir::isolation_options(340), d::SubdivisionPrecision::TrialLocal);
    Observed events;
    cps_stop_reason reason=CPS_STOP_NONE;
    while (reason==CPS_STOP_NONE && run.turns<test_turn_limit) reason=tick(run, events);
    require(reason==CPS_STOP_PRECISION && !run.goal_met, "finite trial cap must give an honest precision stop");
    check_cover(run.cover, roots);
    auto held=run.state();
    require(held.size()==1 && held[0].blocked && held[0].trial_bits>held[0].base_bits,
            "fixture must stop with a raised pending trial separate from base");
    require(events.adopted, "native fixture must exercise successful raised-precision adoption");
    run.engine().reopen(run.options);
    require(run.state()[0].blocked, "same cap must not reopen merely because base is lower");
    run.options.max_work_bits=4096;
    run.engine().reopen(run.options);
    require(!run.state()[0].blocked && run.state()[0].trial_bits==held[0].trial_bits,
            "larger cap must preserve and reopen pending trial");
    while (!run.engine().suspended() && run.turns<test_turn_limit) {
        require(tick(run, events, 0)==CPS_STOP_NONE, "suspension preparation unexpectedly stopped");
        const auto now=run.state();
        require(now.size()==1 && now[0].newton_bits==held[0].newton_bits &&
                    arb_equal(now[0].certificate_radius.value(), held[0].certificate_radius.value()),
                "precision retry must keep e and nominal R fixed");
    }
    require(run.engine().suspended(), "fixture must reach a suspended counting operation");
    held=run.state();
    const auto tests=run.statistics.disk_tests;
    tick(run, events, 0);
    require(run.engine().suspended() && run.state()[0].trial_bits==held[0].trial_bits &&
                run.state()[0].newton_bits==held[0].newton_bits && run.statistics.disk_tests==tests,
            "same trial must resume its counting state without resetting precision/exponent");
    bool cancelled=true;
    run.options.cancel=[](void *flag) { return *static_cast<bool *>(flag) ? 1 : 0; };
    run.options.cancel_data=&cancelled;
    require(tick(run, events)==CPS_STOP_CANCELLED && run.state()[0].trial_bits==held[0].trial_bits,
            "cancellation must preserve pending trial and old cover");
    check_cover(run.cover, roots);
    run.options.cancel=nullptr;
    run.options.cancel_data=nullptr;
    finish(run, events, roots);
    require(events.discarded && events.branched, "native trajectory must exercise post-raise discard and subsequent branching");
    std::cout << "cap/suspend/cancel/resume base=" << held[0].base_bits << " trial=" << held[0].trial_bits
              << " adopted=1 discarded=1 branched=1\n";
}
} // namespace

int main(int argc, char **argv) {
    certpolysolve::ThreadScope lifetime;
    try {
        require(argc==2, "usage: test_subdivision_trial_precision FIXTURE_DIRECTORY");
        const std::filesystem::path directory=argv[1];
        for (const std::string name : {"symmetric-n2-l128", "gaussian-n2-l128"}) full_isolation(directory, name);
        cap_and_suspension(directory);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_subdivision_trial_precision: " << error.what() << '\n';
        return 1;
    }
}
