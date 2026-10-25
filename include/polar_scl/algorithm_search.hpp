#pragma once

#include "polar_scl/portfolio_evaluator.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace polar_scl
{
    struct SelectorSearchSpace
    {
        std::vector<std::size_t> list_size_thresholds{0,2,4};
        std::vector<Rate1Engine> low_list_engines{Rate1Engine::CertifiedFrontier,Rate1Engine::FastSSCLBounded,Rate1Engine::LazyExact};
        std::vector<std::size_t> small_node_thresholds{2,4,8};
        std::vector<double> reliability_gap_thresholds{0.25, 0.5, 1.0, 2.0};
        std::vector<double> parent_metric_gap_thresholds{0.1, 0.25, 0.5, 1.0};
        std::vector<Rate1Engine> small_engines{Rate1Engine::LazyExact,Rate1Engine::FastSSCLBounded,Rate1Engine::CertifiedFrontier};
        std::vector<Rate1Engine> easy_engines{Rate1Engine::CertifiedFrontier,Rate1Engine::LazyExact,Rate1Engine::FastSSCLBounded};
        std::vector<Rate1Engine> hard_engines{Rate1Engine::LazyExact,Rate1Engine::CertifiedFrontier,Rate1Engine::FastSSCLBounded};
        std::vector<std::size_t> frontier_initial_bits{0, 1, 2};
        std::vector<std::size_t> frontier_growth_steps{1, 2};
        std::vector<ExpansionMode> frontier_modes{ExpansionMode::MinimumBound, ExpansionMode::AllCompetitive};
        std::vector<std::size_t> frontier_budgets{64,128,256};
        std::vector<bool> repetition_choices{false,true};
        std::vector<bool> spc_choices{false,true};
    };

    struct SelectorEvolutionConfig
    {
        std::size_t population_size{10};
        std::size_t offspring_per_generation{10};
        std::size_t generations{8};
        std::uint64_t seed{17};
    };

    struct SelectorSearchResult
    {
        std::vector<PortfolioResult> evaluated;
        PortfolioResult best;
        std::vector<double> best_cost_by_generation;
    };

    SelectorSearchResult evolutionary_selector_search(
    const SelectorSearchSpace& space = {},
    const SelectorEvolutionConfig& evolution = {},
    const PortfolioEvaluationConfig& evaluation = {});
}

// namespace polar_scl
