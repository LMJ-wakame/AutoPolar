#pragma once

#include "polar_scl/ir.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace polar_scl
{
    struct EvaluationConfig
    {
        std::vector<std::size_t> node_sizes{8, 16, 32, 64};
        std::vector<std::size_t> list_sizes{2, 4, 8};
        std::size_t trials_per_case{100};
        std::optional<double> quantization_step{};
        std::uint64_t seed{42};
    };

    struct PolicyResult
    {
        Rate1Policy policy;
        std::size_t exact_cases{};
        std::size_t total_cases{};
        double mean_materialized{};
        double p99_materialized{};
        double mean_lazy_materialized{};
        double fallback_rate{};
        double mean_frontier_expansions{};

        bool is_exact() const
        {
            return exact_cases == total_cases;
        }
    };

    PolicyResult evaluate_policy(
    const Rate1Policy& policy,
    const EvaluationConfig& config = {});

    std::vector<PolicyResult> pareto_policies(const std::vector<PolicyResult>& results);
}

// namespace polar_scl
