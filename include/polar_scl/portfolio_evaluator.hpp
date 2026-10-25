#pragma once

#include "polar_scl/algorithm_ir.hpp"
#include "polar_scl/cost.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace polar_scl
{
    struct PortfolioEvaluationConfig
    {
        std::vector<std::size_t> node_sizes{4, 8, 16, 32};
        std::vector<std::size_t> list_sizes{2, 4, 8};
        std::size_t trials_per_case{50};
        std::optional<double> quantization_step{0.25};
        std::uint64_t seed{42};
    };

    struct PortfolioResult
    {
        PolarAlgorithmIR algorithm;
        std::size_t exact_cases{};
        std::size_t total_cases{};
        double mean_proxy_cost{};
        double p99_proxy_cost{};
        double mean_candidates{};
        double fallback_rate{};
        std::vector<std::pair<Rate1Engine, std::size_t>> engine_counts;

        bool is_exact() const
        {
            return exact_cases == total_cases;
        }
    };

    PortfolioResult evaluate_algorithm(
    const PolarAlgorithmIR& algorithm,
    const PortfolioEvaluationConfig& config = {},
    const CostWeights& weights = {});

    PortfolioResult evaluate_selector(
    const SelectorGenome& genome,
    const PortfolioEvaluationConfig& config = {},
    const CostWeights& weights = {});
}

// namespace polar_scl
