#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

namespace polar_scl
{
    enum class ExpansionMode
    {
        MinimumBound,
        AllCompetitive,
    };

    inline std::string to_string(ExpansionMode mode)
    {
        switch(mode)
        {
            case ExpansionMode::MinimumBound: return "minimum_bound";
            case ExpansionMode::AllCompetitive: return "all_competitive";
        }

        throw std::logic_error("unreachable expansion mode");
    }

    struct Rate1Policy
    {
        std::size_t initial_active_bits{1};
        std::size_t growth_step{1};
        ExpansionMode expansion_mode{ExpansionMode::MinimumBound};
        std::size_t candidate_budget{4096};

        auto operator<=>(const Rate1Policy&) const = default;
    };

    inline void validate_policy(const Rate1Policy& policy)
    {
        if(policy.expansion_mode != ExpansionMode::MinimumBound && policy.expansion_mode != ExpansionMode::AllCompetitive)
            throw std::invalid_argument("unsupported expansion mode");
        if(policy.growth_step < 1)
        {
            throw std::invalid_argument("growth_step must be at least one");
        }

        if(policy.candidate_budget < 1)
        {
            throw std::invalid_argument("candidate_budget must be positive");
        }
    }
}

// namespace polar_scl
