#pragma once
// The program route's second version (docs/slp-second-version-note-20261017.md): the roots of a program certified
// from the program's own evaluations in balls, without the regeneration's n^2 node products. Phase 2a: separated
// roots by Newton inclusion disks in a disjoint cover. Phase 2b: clusters by root counts on circles (the argument
// principle), and a root bound from such counts.
//
// The Newton disk. For a polynomial p of degree n, p'(z) / p(z) is the sum of 1 / (z - root) over the n roots
// counted with multiplicity, so some root lies within n |p(z) / p'(z)| of z; in balls the radius is n times an
// upper bound of |p(z)| over a lower bound of |p'(z)|, both evaluated at the exact point z. n pairwise disjoint
// disks of a polynomial of degree n, each holding at least one root, hold exactly one each, and every root is
// simple (a multiple root would leave fewer than n distinct roots for n disjoint disks). So the n disks, published
// as squares and found pairwise disjoint, are a complete cover: count 1 and multiplicity 1 in every region.
//
// Only the whole cover certifies. A disk disjoint from all the others proves nothing alone: when two points have
// converged to one root, another root has no point near it and may lie in some disjoint disk, which then holds two
// roots. So the pass certifies when all n disks are disjoint and meet the goal.
//
// Counted clusters (phase 2b). When disks keep overlapping (a multiple root, or roots closer than the goal), the
// overlapping points are grouped into components, and each component's circle gets a certified root count m. Its
// region, a square holding the circle, then holds at least m roots, and each separated disk at least one; regions
// that are pairwise disjoint hold distinct roots, so when the counts and the separated disks add up to n, every
// region holds exactly its count and no root lies outside them. Otherwise the pass fails, and the first version's
// refinement is the fallback.
#include "core/cover.hpp"
#include "core/query.hpp"
#include "core/engine_statistics.hpp"
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace certpolysolve::detail {
class MonicSource;

class SeparatedRootsPass {
  public:
    enum class State : uint8_t { Running, Certified, Failed };
    // From the points' exact midpoints, one a root of a source of degree n, at `precision`.
    SeparatedRootsPass(acb_srcptr points, size_t n, slong precision);
    // One round: each point whose disk is not current evaluated once (p and p' in balls at the pass's
    // precision) and its disk formed; the overlap sweep over all n regions; an Aberth step for each point whose
    // disk overlaps a disk no larger than its own, a Newton step for each that is only too wide for the goal,
    // unless its value is at its rounding (then the precision doubles, to the cap). Certified when every disk is
    // current, the n regions are pairwise disjoint and each meets the goal. When the goal admits counted clusters,
    // the cluster step (cluster_cover) runs once when every region meets the goal and only overlaps remain, at the
    // first round without progress, before re-homing starts (two
    // points at a multiple root stall there, and no free root exists for re-homing to find), and again where the
    // pass would fail with overlapping disks (after the rounds allowed, when the overlaps stop shrinking, or at the
    // cap). A cancelled early attempt is not repeated: the stall count has moved on, and the next round re-homes.
    // Failed when a derivative has no correct bits at the cap, or when the last cluster step does not certify.
    // Cancellation leaves the state Running, the points as they were, and the round to do again.
    // `undivided`, when given, is z^v times the source's polynomial for some v >= 1: a program's polynomial p before
    // its root of order v at zero was divided out (the source is q = p / z^v). Off zero it has q's roots, and its
    // enclosures over balls are tighter: q' = (p' - v p / z) / z^v carries p's value, and a ball's p loses the
    // cancellation of the program's last sums ((P + Q) - Q), which makes q'(B) far wider than p'(B). The cluster step
    // counts on it in place of the source on a circle that lies far from zero against its radius
    // (heuristics::program_undivided_count_margin), where its count is q's, and on the source otherwise.
    // History: docs/slp-second-version-note-20261017.md, section 10.
    State advance(MonicSource &source, const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics,
                  MonicSource *undivided = nullptr);
    State state() const { return state_; }
    // The cover, meaningful when Certified: one region a point (count and multiplicity 1), or the separated
    // points' disks and the counted clusters (count m; multiplicity 1 when m is 1, unknown otherwise).
    const std::vector<Root> &cover() const { return cover_; }
    uint64_t rounds() const { return rounds_; }
    uint64_t evaluations() const { return evaluations_; }
    uint64_t raises() const { return raises_; }
    uint64_t steps() const { return steps_; }
    // Points moved off a root another point held (re-homed).
    uint64_t rehomed() const { return rehomed_; }
    // Counted clusters in the certified cover, and the root counts on circles taken (each with its arcs).
    uint64_t clusters() const { return clusters_; }
    uint64_t counts() const { return counts_; }
    uint64_t count_arcs() const { return count_arcs_; }
    slong precision() const { return precision_; }
    // After the last round: points whose disk overlaps a disk no larger than its own (the culprits), regions too
    // wide for the goal.
    size_t overlapping() const { return overlapping_; }
    size_t wide() const { return wide_; }
    const std::string &failure() const { return failure_; }

  private:
    enum class ClusterOutcome : uint8_t { Certified, None, Cancelled };
    // The cluster step: the components of the overlapping pairs, a certified count on a circle about each (the
    // precision doubling once when a circle's values are unreadable), each circle shrunk to the goal about its
    // points' centroid (the precision doubling, to the cap, while the goal's circle is unreadable; a readable
    // count that differs ends the step: distinct roots, not one cluster), then the cover's completeness,
    // disjointness and goal checked. A count on `undivided` where advance() says.
    ClusterOutcome cluster_cover(MonicSource &source, MonicSource *undivided, const cps_options &options,
                                 certpolysolve::detail::CandidateStatistics &statistics,
                                 const std::vector<std::pair<size_t, size_t>> &pairs,
                                 const std::vector<uint8_t> &involved);
    std::vector<ComplexBall> points_, corrections_;
    std::vector<Root> regions_, cover_;
    std::vector<uint8_t> current_, limited_, unreadable_;
    State state_ = State::Running;
    slong precision_ = 0;
    uint64_t rounds_ = 0, evaluations_ = 0, raises_ = 0, steps_ = 0, stalled_rounds_ = 0, rehomed_ = 0;
    uint64_t clusters_ = 0, counts_ = 0, count_arcs_ = 0;
    size_t overlapping_ = 0, wide_ = 0, best_overlapping_ = SIZE_MAX;
    double best_widest_ = std::numeric_limits<double>::infinity(); // log2, the widest culprit disk's best so far
    bool raised_since_progress_ = false; // a raise made and no progress since: no further raise until progress
    bool converged_tried_ = false;       // the cluster step tried once every region met the goal
    std::string failure_;
};

} // namespace certpolysolve::detail
