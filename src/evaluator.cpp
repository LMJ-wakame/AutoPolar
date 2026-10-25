#include "polar_scl/evaluator.hpp"
#include "polar_scl/rate1.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>

namespace polar_scl
{
    namespace
    {
        double quantize(double value, const std::optional<double>& step)
        {
            if(!step) return value;
            return std::nearbyint(value / *step) * *step;
        }

        double mean(const std::vector<double>& values)
        {
            return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
        }

        double p99(std::vector<double> values)
        {
            if(values.empty()) return 0.0;
            std::sort(values.begin(), values.end());
            std::size_t idx = static_cast<std::size_t>(0.99 * static_cast<double>(values.size()));
            if(idx > 0) --idx;
            idx = std::min(idx, values.size() - 1);
            return values[idx];
        }

        bool same_candidates(const std::vector<Rate1Candidate>& a, const std::vector<Rate1Candidate>& b)
        {
            if(a.size() != b.size()) return false;
            for(std::size_t i = 0; i < a.size(); ++i) if(!(local_rank_key(a[i]) == local_rank_key(b[i]))) return false;
            return true;
        }
    }

    // namespace

    PolicyResult evaluate_policy(const Rate1Policy& policy, const EvaluationConfig& cfg)
    {
        std::mt19937_64 rng(cfg.seed);
        std::exponential_distribution<double> metric_dist(1.0);
        std::normal_distribution<double> llr_dist(0.0, 3.0);
        std::size_t exact = 0, fallbacks = 0, total = 0;
        std::vector<double> mats, lazy_mats, expansions;
        for(auto node_size : cfg.node_sizes) for(auto list_size : cfg.list_sizes) for(std::size_t t = 0; t < cfg.trials_per_case; ++t)
        {
            std::vector<ParentPath> parents;
            for(std::size_t id = 0; id < list_size; ++id)
            {
                std::vector<double> llrs(node_size);
                for(double& v : llrs) v = quantize(llr_dist(rng), cfg.quantization_step);
                parents.push_back({id, quantize(metric_dist(rng), cfg.quantization_step), std::move(llrs)});
            }

            auto expected = lazy_exact_rate1_top_l(parents, list_size);
            auto actual = adaptive_rate1_top_l(parents, list_size, policy);
            exact += same_candidates(expected.candidates, actual.candidates) ? 1 : 0;
            mats.push_back(static_cast<double>(actual.stats.candidates_materialized));
            lazy_mats.push_back(static_cast<double>(expected.stats.candidates_materialized));
            expansions.push_back(static_cast<double>(actual.stats.frontier_expansions));
            fallbacks += actual.stats.fallback_used ? 1 : 0;
            ++total;
        }

        return {policy, exact, total, mean(mats), p99(mats), mean(lazy_mats), static_cast<double>(fallbacks)/static_cast<double>(total), mean(expansions)};
    }

    std::vector<PolicyResult> pareto_policies(const std::vector<PolicyResult>& results)
    {
        std::vector<PolicyResult> frontier;
        for(const auto& c : results)
        {
            if(!c.is_exact()) continue;
            bool dominated = false;
            for(const auto& o : results)
            {
                if(!o.is_exact()) continue;
                if(o.mean_materialized <= c.mean_materialized && o.fallback_rate <= c.fallback_rate &&
                (o.mean_materialized < c.mean_materialized || o.fallback_rate < c.fallback_rate))
                {
                    dominated = true;
                    break;
                }
            }

            if(!dominated) frontier.push_back(c);
        }

        std::sort(frontier.begin(), frontier.end(), [](const auto& a, const auto& b)
        {
            return std::tie(a.mean_materialized, a.fallback_rate) < std::tie(b.mean_materialized, b.fallback_rate);
        }

        );
        return frontier;
    }
}

// namespace polar_scl
