#include "polar_scl/selector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace polar_scl
{
    namespace
    {
        double mean(const std::vector<double>& values)
        {
            if(values.empty()) return 0.0;
            return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
        }

        double median_sorted(const std::vector<double>& values)
        {
            const std::size_t n=values.size();
            if(n%2==1)
            {
                return values[n/2];
            }
            return (values[n/2-1]+values[n/2])/2.0;
        }
    }

    // namespace

    Rate1NodeState build_rate1_state(const std::vector<ParentPath>& parents, std::size_t list_size)
    {
        if(parents.empty()) throw std::invalid_argument("at least one parent is required");
        const auto node_size = parents.front().llrs.size();
        std::vector<double> firsts, seconds, medians, gaps;
        for(const auto& parent : parents)
        {
            std::vector<double> ordered;
            ordered.reserve(parent.llrs.size());
            for(double v : parent.llrs) ordered.push_back(std::abs(v));
            std::sort(ordered.begin(), ordered.end());
            firsts.push_back(ordered.front());
            seconds.push_back(ordered.size() > 1 ? ordered[1] : ordered.front());
            medians.push_back(median_sorted(ordered));
            const std::size_t pivot = std::min(std::max(list_size > 0 ? list_size - 1 : 0, std::size_t{1}), ordered.size() - 1);
            gaps.push_back(std::max(0.0, ordered[pivot] - ordered[pivot - 1]));
        }

        std::vector<double> metrics;
        for(const auto& p : parents) metrics.push_back(p.path_metric);
        std::sort(metrics.begin(), metrics.end());
        const double metric_gap = metrics.size() > 1 ? metrics[1] - metrics[0] : std::numeric_limits<double>::infinity();
        return {node_size, list_size, parents.size(), mean(firsts), mean(seconds), mean(gaps), metric_gap, mean(medians)};
    }

    Rate1NodeIR select_rate1_ir(const SelectorGenome& genome, const Rate1NodeState& state)
    {
        validate_selector_genome(genome);
        if(genome.list_size_threshold>0&&state.list_size<=genome.list_size_threshold)
        {
            return genome.ir_for_engine(genome.low_list_engine);
        }

        if(state.node_size<=genome.small_node_threshold)
        {
            return genome.ir_for_engine(genome.small_engine);
        }

        const bool easy=state.reliability_gap>=genome.reliability_gap_threshold||
            state.parent_metric_gap>=genome.parent_metric_gap_threshold;
        return genome.ir_for_engine(easy?genome.easy_engine:genome.hard_engine);
    }
}

// namespace polar_scl
