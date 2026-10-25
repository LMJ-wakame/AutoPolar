#pragma once

#include "polar_scl/evaluator.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace polar_scl
{
    struct SearchSpace
    {
        std::vector<std::size_t> initial_active_bits{0, 1, 2, 3};
        std::vector<std::size_t> growth_steps{1, 2, 4};
        std::vector<ExpansionMode> expansion_modes{ExpansionMode::MinimumBound, ExpansionMode::AllCompetitive};
        std::vector<std::size_t> candidate_budgets{32, 64, 128, 256};

        std::vector<Rate1Policy> enumerate_policies() const;
        void validate() const;
    };

    struct EvolutionConfig
    {
        std::size_t population_size{8};
        std::size_t offspring_per_generation{8};
        std::size_t generations{6};
        std::uint64_t seed{7};
    };

    struct SearchResult
    {
        std::vector<PolicyResult> evaluated;
        std::vector<PolicyResult> pareto_front;
        std::vector<double> best_mean_by_generation;
    };

    SearchResult evolutionary_search(
    const SearchSpace& space = {},
    const EvolutionConfig& evolution = {},
    const EvaluationConfig& evaluation = {});
}

// namespace polar_scl
