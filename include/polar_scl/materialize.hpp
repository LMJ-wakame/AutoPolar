#pragma once

#include "polar_scl/algorithm_ir.hpp"
#include "polar_scl/rate1.hpp"

#include <optional>

namespace polar_scl
{
    Rate1SearchResult execute_rate1_ir(
    const Rate1NodeIR& ir,
    const std::vector<ParentPath>& parents,
    std::size_t list_size);

    struct IntegratedRate1Result
    {
        std::optional<std::vector<Rate1Candidate>> candidates;
        SearchStats stats;
    };

    IntegratedRate1Result execute_rate1_integrated(
    const Rate1NodeIR& ir,
    const std::vector<ParentPath>& parents,
    std::size_t list_size,
    CandidateRankFn rank_key,
    HiddenRankFn hidden_rank_key,
    bool enable_cutoff_tie_guard=true);
}

// namespace polar_scl
