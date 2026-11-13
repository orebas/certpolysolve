// Factor statistics aggregation: ordinary unsigned counters, maxima and counted
// work retain their existing arithmetic. Local axis/prediction diagnostics and
// the current round's shortfall histogram stay with their state. Called by the
// prepared factor loop (prepared_piecewise.cpp); no scheduling or publication here.
#include "methods/piecewise/piecewise_statistics.hpp"
#include "core/statistics_accumulation.hpp"
#include <algorithm>

namespace certpolysolve::detail {
PiecewiseRootStatistics &PiecewiseRootStatistics::operator+=(const PiecewiseRootStatistics &other) {
    work += other.work;
    finishing_work += other.finishing_work;
    statistics::add_ordinary(certificate, other.certificate, statistics::chart_sum_fields);
    statistics::add_ordinary(certificate.model_error, other.certificate.model_error, statistics::model_error_sum_fields);
    statistics::add_ordinary(*this, other, statistics::piecewise_sum_fields);
    for (const auto field : statistics::piecewise_max_fields)
        this->*field = std::max(this->*field, other.*field);
    return *this;
}
} // namespace certpolysolve::detail
