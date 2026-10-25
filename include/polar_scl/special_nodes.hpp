#pragma once

#include "polar_scl/rate1.hpp"

#include <cstddef>
#include <vector>

namespace polar_scl
{
    Rate1SearchResult repetition_node_top_l(
        const std::vector<ParentPath>& parents,
        std::size_t list_size,
        CandidateRankFn rank_key={});

    Rate1SearchResult spc_node_top_l(
        const std::vector<ParentPath>& parents,
        std::size_t list_size,
        CandidateRankFn rank_key={});
}
