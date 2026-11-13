#pragma once
// The synchronous query envelope: validation, defaults and cancellation latching.
// It owns the callback's state on the stack; moving/copying it would invalidate
// next.cancel_data, so both are forbidden. Routes borrow next for synchronous execution.
#include "solve/query_input.hpp"
#include <limits>
#include <optional>
#include <utility>

namespace certpolysolve::detail {
class QueryExecution {
    struct Cancellation {
        CancelCallback callback = nullptr;
        void *data = nullptr;
        bool observed = false;
        static int poll(void *pointer) {
            auto &self = *static_cast<Cancellation *>(pointer);
            if (!self.observed && self.callback && self.callback(self.data)) self.observed = true;
            return self.observed;
        }
    } cancellation_;
  public:
    cps_options next;
    uint64_t degree = 0;
    bool zero_input = false;
    std::unique_ptr<RealRegionFilter> region;
    std::optional<cps_result> refusal;
    QueryExecution(const QueryInput &source, const cps_options *proposed);
    QueryExecution(const QueryExecution &) = delete;
    QueryExecution &operator=(const QueryExecution &) = delete;
    QueryExecution(QueryExecution &&) = delete;
    QueryExecution &operator=(QueryExecution &&) = delete;
    template <class Route> cps_result run(Route &&route) {
        try {
            return std::forward<Route>(route)();
        } catch (const std::bad_alloc &) {
            throw;
        } catch (const ExpansionRefused &) {
            throw;
        } catch (const std::exception &error) {
            throw std::runtime_error(std::string("numerical solve failed: ") + error.what());
        }
    }
};
} // namespace certpolysolve::detail
