#include "polar_scl/cost.hpp"

namespace polar_scl
{
    double ComplexityVector::proxy_cost(const CostWeights& w) const
    {
        return w.candidate * static_cast<double>(candidates)
        + w.heap_operation * static_cast<double>(heap_operations())
        + w.frontier_expansion * static_cast<double>(frontier_expansions)
        + w.certificate_check * static_cast<double>(certificate_checks);
    }

    ComplexityVector complexity_from_stats(const SearchStats& s)
    {
        return {s.candidates_materialized,
            s.subset_states_pushed + s.global_heap_pushes,
            s.subset_states_popped + s.global_heap_pops,
            s.frontier_expansions,
            s.certificate_checks,
            s.fallback_used};
    }
}

// namespace polar_scl
