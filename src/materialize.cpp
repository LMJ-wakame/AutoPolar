#include "polar_scl/materialize.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace polar_scl
{
    namespace
    {
        bool metric_close(double lhs,double rhs,double tolerance=1e-9)
        {
            return std::abs(lhs-rhs)<=tolerance*std::max({1.0,std::abs(lhs),std::abs(rhs)});
        }
    }

    Rate1SearchResult execute_rate1_ir(const Rate1NodeIR& ir, const std::vector<ParentPath>& parents, std::size_t list_size)
    {
        validate_rate1_ir(ir);
        switch(ir.engine)
        {
            case Rate1Engine::Exhaustive: return exhaustive_rate1_top_l(parents, list_size);
            case Rate1Engine::LazyExact: return lazy_exact_rate1_top_l(parents, list_size);
            case Rate1Engine::FastSSCLBounded: return fast_sscl_bounded_rate1_top_l(parents, list_size);
            case Rate1Engine::CertifiedFrontier: return adaptive_rate1_top_l(parents, list_size, ir.frontier_policy());
        }

        throw std::logic_error("unreachable Rate-1 engine");
    }

    IntegratedRate1Result execute_rate1_integrated(
    const Rate1NodeIR& ir,
    const std::vector<ParentPath>& parents,
    std::size_t list_size,
    CandidateRankFn rank_key,
    HiddenRankFn hidden_rank_key,
    bool enable_cutoff_tie_guard)
    {
        validate_rate1_ir(ir);
        if(ir.engine == Rate1Engine::CertifiedFrontier)
        {
            auto result = adaptive_rate1_top_l(parents, list_size, ir.frontier_policy(), std::move(rank_key), std::move(hidden_rank_key));
            return {std::move(result.candidates), std::move(result.stats)};
        }

        const std::size_t probe_size = list_size + 1;
        Rate1SearchResult probe;
        switch(ir.engine)
        {
            case Rate1Engine::Exhaustive: probe = exhaustive_rate1_top_l(parents, probe_size);
            break;
            case Rate1Engine::LazyExact: probe = lazy_exact_rate1_top_l(parents, probe_size);
            break;
            case Rate1Engine::FastSSCLBounded: probe = fast_sscl_bounded_rate1_top_l(parents, probe_size);
            break;
            case Rate1Engine::CertifiedFrontier: break;
        }

        if(enable_cutoff_tie_guard&&probe.candidates.size()>list_size&&
        metric_close(
            probe.candidates[list_size-1].path_metric,
            probe.candidates[list_size].path_metric))
        {
            probe.stats.fallback_used = true;
            probe.stats.fallback_reason = "whole-frame tie at Top-L cutoff";
            return {std::nullopt, std::move(probe.stats)};
        }

        if(probe.candidates.size() > list_size) probe.candidates.resize(list_size);
        std::sort(probe.candidates.begin(), probe.candidates.end(), [&](const auto& a, const auto& b)
        {
            return rank_key(a) < rank_key(b);
        }

        );
        return {std::move(probe.candidates), std::move(probe.stats)};
    }
}

// namespace polar_scl
