#pragma once

#include "polar_scl/algorithm_ir.hpp"
#include "polar_scl/rate1.hpp"

#include <cstddef>
#include <vector>

namespace polar_scl
{
    struct Rate1NodeState
    {
        std::size_t node_size{};
        std::size_t list_size{};
        std::size_t parent_count{};
        double min_reliability{};
        double second_reliability{};
        double reliability_gap{};
        double parent_metric_gap{};
        double median_reliability{};
    };

    Rate1NodeState build_rate1_state(
    const std::vector<ParentPath>& parents,
    std::size_t list_size);

    Rate1NodeIR select_rate1_ir(
    const SelectorGenome& genome,
    const Rate1NodeState& state);
}

// namespace polar_scl
