#include "polar_scl/portfolio_evaluator.hpp"
#include "polar_scl/portfolio.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <random>

namespace polar_scl
{
    namespace
    {
        double quantize(double value,const std::optional<double>&step)
        {
            return step?std::round(value/(*step))*(*step):value;
        }

        double mean(const std::vector<double>&values)
        {
            if(values.empty())
            {
                return 0.0;
            }
            return std::accumulate(values.begin(),values.end(),0.0)/static_cast<double>(values.size());
        }

        double p99(std::vector<double>values)
        {
            if(values.empty())
            {
                return 0.0;
            }
            std::sort(values.begin(),values.end());
            std::size_t index=static_cast<std::size_t>(0.99*static_cast<double>(values.size()));
            if(index>0)
            {
                --index;
            }
            return values[std::min(index,values.size()-1)];
        }

        bool same(const std::vector<Rate1Candidate>&lhs,const std::vector<Rate1Candidate>&rhs)
        {
            if(lhs.size()!=rhs.size())
            {
                return false;
            }
            for(std::size_t index=0;index<lhs.size();++index)
            {
                if(!(local_rank_key(lhs[index])==local_rank_key(rhs[index])))
                {
                    return false;
                }
            }
            return true;
        }
    }

    PortfolioResult evaluate_algorithm(const PolarAlgorithmIR&algorithm,const PortfolioEvaluationConfig&config,const CostWeights&weights)
    {
        std::mt19937_64 rng(config.seed);
        std::exponential_distribution<double>metric_distribution(1.0);
        std::normal_distribution<double>llr_distribution(0.0,3.0);
        std::vector<double>costs;
        std::vector<double>candidates;
        std::size_t exact=0;
        std::size_t fallbacks=0;
        std::size_t total=0;
        std::map<Rate1Engine,std::size_t>engines;

        for(auto node_size:config.node_sizes)
        {
            for(auto list_size:config.list_sizes)
            {
                for(std::size_t trial=0;trial<config.trials_per_case;++trial)
                {
                    std::vector<ParentPath>parents;
                    for(std::size_t parent_id=0;parent_id<list_size;++parent_id)
                    {
                        std::vector<double>llrs(node_size);
                        for(double&llr:llrs)
                        {
                            llr=quantize(llr_distribution(rng),config.quantization_step);
                        }
                        parents.push_back({parent_id,quantize(metric_distribution(rng),config.quantization_step),std::move(llrs)});
                    }

                    auto expected=lazy_exact_rate1_top_l(parents,list_size);
                    auto actual=portfolio_rate1_top_l(algorithm,parents,list_size);
                    exact+=same(actual.candidates,expected.candidates)?1:0;

                    auto complexity=complexity_from_stats(actual.stats);
                    costs.push_back(complexity.proxy_cost(weights));
                    candidates.push_back(static_cast<double>(complexity.candidates));
                    fallbacks+=complexity.fallback_used?1:0;
                    ++engines[actual.ir.engine];
                    ++total;
                }
            }
        }

        std::vector<std::pair<Rate1Engine,std::size_t>>engine_counts(engines.begin(),engines.end());
        return {
            algorithm,
            exact,
            total,
            mean(costs),
            p99(costs),
            mean(candidates),
            total==0?0.0:static_cast<double>(fallbacks)/static_cast<double>(total),
            std::move(engine_counts)
        };
    }

    PortfolioResult evaluate_selector(const SelectorGenome&genome,const PortfolioEvaluationConfig&config,const CostWeights&weights)
    {
        PolarAlgorithmIR ir;
        ir.selector_mode=SelectorMode::Rule;
        ir.selector=genome;
        return evaluate_algorithm(ir,config,weights);
    }
}
