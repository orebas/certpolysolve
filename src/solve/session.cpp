// Public call dispatch. Validation/defaults/cancellation belong to QueryExecution;
// the real-only route and ordinary engine orchestration are independent siblings.
#include "solve/solve_state.hpp"
#include "solve/query_execution.hpp"
#include "solve/complex_query.hpp"
#include "solve/real_only_route.hpp"

namespace certpolysolve::detail {
struct SessionState::Impl {
    QueryInput source;
    ComplexQueryState complex;
    std::unique_ptr<RealOnlyRoute> real_route;
    bool advanced = false;
    explicit Impl(const DensePolynomial &p) : source(p), complex(source) {}
    explicit Impl(DensePolynomial &&p) : source(std::move(p)), complex(source) {}
    explicit Impl(const cps_polynomial &p) : source(p), complex(source) {}
    explicit Impl(const cps_secular &p) : source(p), complex(source) {}
    explicit Impl(const cps_program &p) : source(p), complex(source) {}
    cps_result advance(const cps_options *proposed) {
        if (advanced) throw std::logic_error("a solve state is advanced once");
        advanced = true;
        QueryExecution query(source, proposed);
        if (query.refusal) return std::move(*query.refusal);
        return query.run([&] {
            const auto &next = query.next;
            if (source.polynomial() && source.real_route_applies(next) &&
                (next.goal == CPS_ROOT_ACCURACY || next.goal == CPS_ROOT_APPROXIMATION ||
                 next.goal == CPS_ROOT_ISOLATION)) {
                real_route = source.input
                    ? std::make_unique<RealOnlyRoute>(std::move(*source.input), next, query.degree,
                                                       query.zero_input, std::move(query.region))
                    : std::make_unique<RealOnlyRoute>(*source.terms, next, query.degree,
                                                       query.zero_input, std::move(query.region));
                source.input.reset();
                return real_route->advance(next);
            }
            return complex.advance(query.next, query.degree);
        });
    }
};
SessionState::SessionState(const DensePolynomial &p) : impl_(std::make_unique<Impl>(p)) {}
SessionState::SessionState(DensePolynomial &&p) : impl_(std::make_unique<Impl>(std::move(p))) {}
SessionState::SessionState(const cps_polynomial &borrowed) : impl_(std::make_unique<Impl>(borrowed)) {}
SessionState::SessionState(const cps_secular &p) : impl_(std::make_unique<Impl>(p)) {}
SessionState::SessionState(const cps_program &p) : impl_(std::make_unique<Impl>(p)) {}
SessionState::~SessionState() = default;
cps_result SessionState::advance(const cps_options *options) { return impl_->advance(options); }
} // namespace certpolysolve::detail
