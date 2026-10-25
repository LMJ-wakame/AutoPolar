#pragma once

#include "polar_scl/ir.hpp"

#include <cstddef>
#include <string>

namespace polar_scl
{
    enum class Rate1Engine
    {
        Exhaustive,
        LazyExact,
        FastSSCLBounded,
        CertifiedFrontier,
    };

    enum class SelectorMode
    {
        Fixed,
        Rule,
    };

    std::string to_string(Rate1Engine engine);

    struct Rate1NodeIR
    {
        Rate1Engine engine{Rate1Engine::CertifiedFrontier};
        std::size_t initial_active_bits{1};
        std::size_t growth_step{1};
        ExpansionMode expansion_mode{ExpansionMode::MinimumBound};
        std::size_t candidate_budget{256};

        Rate1Policy frontier_policy() const;
        auto operator<=>(const Rate1NodeIR&) const = default;
    };

    struct SelectorGenome
    {
        std::size_t list_size_threshold{4};
        Rate1Engine low_list_engine{Rate1Engine::CertifiedFrontier};
        std::size_t small_node_threshold{4};
        double reliability_gap_threshold{1.0};
        double parent_metric_gap_threshold{0.5};
        Rate1Engine small_engine{Rate1Engine::Exhaustive};
        Rate1Engine easy_engine{Rate1Engine::CertifiedFrontier};
        Rate1Engine hard_engine{Rate1Engine::LazyExact};
        std::size_t frontier_initial_active_bits{0};
        std::size_t frontier_growth_step{1};
        ExpansionMode frontier_expansion_mode{ExpansionMode::MinimumBound};
        std::size_t frontier_candidate_budget{256};
        bool enable_repetition{true};
        bool enable_spc{true};

        Rate1NodeIR ir_for_engine(Rate1Engine engine) const;
        auto operator<=>(const SelectorGenome&) const = default;
    };

    struct PolarAlgorithmIR
    {
        SelectorMode selector_mode{SelectorMode::Rule};
        Rate1NodeIR fixed_rate1{};
        SelectorGenome selector{};
    };

    void validate_rate1_ir(const Rate1NodeIR& ir);
    void validate_selector_genome(const SelectorGenome& genome);
    void validate_algorithm_ir(const PolarAlgorithmIR& ir);
}

// namespace polar_scl
