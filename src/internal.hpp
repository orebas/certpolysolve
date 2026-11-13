#pragma once
#include <certpolysolve/certpolysolve.h>
#include "core/engine_statistics.hpp"
#include "core/number.hpp"
#include "core/heuristics.hpp"
#include "core/engine_config.hpp"
#include "core/resource_limits.hpp"
#include "numeric/exact/prepared_structure.hpp"
#include <array>
#include <atomic>
#include <chrono>
#include <map>
#include <stdexcept>
#include <memory>
#include <vector>
#include <optional>
#include <tuple>
#include <functional>
// The solver's internal declarations as one include, in five parts: the cover (cover.hpp), the polynomial as a
// solve receives it (polynomial_input.hpp), a call's options and result (query.hpp), the engines' switches
// (engine_policy.hpp), and the solve's entry points and states (solve_state.hpp). The library's own files include the
// parts they use; the tests and the benchmark tools include this file.
#include "core/cover.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "core/engine_policy.hpp"
#include "solve/solve_state.hpp"
