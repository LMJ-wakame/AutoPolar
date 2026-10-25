#include "polar_scl/algorithm_search.hpp"

#include <algorithm>
#include <map>
#include <random>
#include <stdexcept>
#include <tuple>

namespace polar_scl
{
    namespace
    {
        template<class T> const T& choose(const std::vector<T>& v,std::mt19937_64&r)
        {
            std::uniform_int_distribution<std::size_t>d(0,v.size()-1);
            return v[d(r)];
        }

        SelectorGenome random_selector(const SelectorSearchSpace&s,std::mt19937_64&r)
        {
            SelectorGenome g;
            g.list_size_threshold=choose(s.list_size_thresholds,r);
            g.low_list_engine=choose(s.low_list_engines,r);
            g.small_node_threshold=choose(s.small_node_thresholds,r);
            g.reliability_gap_threshold=choose(s.reliability_gap_thresholds,r);
            g.parent_metric_gap_threshold=choose(s.parent_metric_gap_thresholds,r);
            g.small_engine=choose(s.small_engines,r);
            g.easy_engine=choose(s.easy_engines,r);
            g.hard_engine=choose(s.hard_engines,r);
            g.frontier_initial_active_bits=choose(s.frontier_initial_bits,r);
            g.frontier_growth_step=choose(s.frontier_growth_steps,r);
            g.frontier_expansion_mode=choose(s.frontier_modes,r);
            g.frontier_candidate_budget=choose(s.frontier_budgets,r);
            validate_selector_genome(g);
            return g;
        }

        SelectorGenome mutate(const SelectorGenome&g,const SelectorSearchSpace&s,std::mt19937_64&r)
        {
            std::uniform_int_distribution<int> sd(0,11);
            int slot=sd(r);
            auto c=g;
            auto diff=[&](const auto& vals,const auto& cur)
            {
                using T=std::decay_t<decltype(cur)>;
                std::vector<T>x;
                for(auto v:vals)if(v!=cur)x.push_back(v);
                return x.empty()?cur:choose(x,r);
            };
            switch(slot)
            {
                case 0:c.list_size_threshold=diff(s.list_size_thresholds,g.list_size_threshold);
                break;
                case 1:c.low_list_engine=diff(s.low_list_engines,g.low_list_engine);
                break;
                case 2:c.small_node_threshold=diff(s.small_node_thresholds,g.small_node_threshold);
                break;
                case 3:c.reliability_gap_threshold=diff(s.reliability_gap_thresholds,g.reliability_gap_threshold);
                break;
                case 4:c.parent_metric_gap_threshold=diff(s.parent_metric_gap_thresholds,g.parent_metric_gap_threshold);
                break;
                case 5:c.small_engine=diff(s.small_engines,g.small_engine);
                break;
                case 6:c.easy_engine=diff(s.easy_engines,g.easy_engine);
                break;
                case 7:c.hard_engine=diff(s.hard_engines,g.hard_engine);
                break;
                case 8:c.frontier_initial_active_bits=diff(s.frontier_initial_bits,g.frontier_initial_active_bits);
                break;
                case 9:c.frontier_growth_step=diff(s.frontier_growth_steps,g.frontier_growth_step);
                break;
                case 10:c.frontier_expansion_mode=diff(s.frontier_modes,g.frontier_expansion_mode);
                break;
                case 11:c.frontier_candidate_budget=diff(s.frontier_budgets,g.frontier_candidate_budget);
                break;
            }

            validate_selector_genome(c);
            return c;
        }

        auto key(const PortfolioResult&r)
        {
            return std::tuple{!r.is_exact(),r.mean_proxy_cost,r.p99_proxy_cost,r.mean_candidates,r.fallback_rate};
        }
    }

    // namespace

    SelectorSearchResult evolutionary_selector_search(const SelectorSearchSpace&space,const SelectorEvolutionConfig&evo,const PortfolioEvaluationConfig&evaluation)
    {
        if(evo.population_size<1||evo.offspring_per_generation<1||evo.generations<1)throw std::invalid_argument("population/offspring/generations must be positive");
        std::mt19937_64 rng(evo.seed);
        std::map<SelectorGenome,PortfolioResult> cache;
        auto score=[&](const SelectorGenome&g)->const PortfolioResult&
        {
            auto[it,ins]=cache.try_emplace(g);
            if(ins)it->second=evaluate_selector(g,evaluation);
            return it->second;
        };
        std::vector<SelectorGenome> seeds;
        {
            SelectorGenome g;
            g.list_size_threshold=0;
            g.small_node_threshold=2;
            g.small_engine=g.easy_engine=g.hard_engine=Rate1Engine::LazyExact;
            seeds.push_back(g);
        }

        {
            SelectorGenome g;
            g.list_size_threshold=0;
            g.small_node_threshold=2;
            g.reliability_gap_threshold=.25;
            g.parent_metric_gap_threshold=.1;
            g.small_engine=Rate1Engine::LazyExact;
            g.easy_engine=g.hard_engine=Rate1Engine::CertifiedFrontier;
            g.frontier_initial_active_bits=0;
            g.frontier_candidate_budget=256;
            seeds.push_back(g);
        }

        {
            SelectorGenome g;
            g.list_size_threshold=0;
            g.small_node_threshold=4;
            g.small_engine=Rate1Engine::LazyExact;
            g.easy_engine=Rate1Engine::CertifiedFrontier;
            g.hard_engine=Rate1Engine::LazyExact;
            g.frontier_initial_active_bits=0;
            seeds.push_back(g);
        }

        std::vector<SelectorGenome> pop;
        for(auto g:seeds)if(pop.size()<evo.population_size&&std::find(pop.begin(),pop.end(),g)==pop.end())pop.push_back(g);
        while(pop.size()<evo.population_size)
        {
            auto g=random_selector(space,rng);
            if(std::find(pop.begin(),pop.end(),g)==pop.end())pop.push_back(g);
        }

        std::vector<double> history;
        for(std::size_t gen=0;gen<evo.generations;++gen)
        {
            std::vector<PortfolioResult> parents;
            for(auto&g:pop)parents.push_back(score(g));
            std::sort(parents.begin(),parents.end(),[](auto&a,auto&b)
            {
                return key(a)<key(b);
            }

            );
            std::vector<SelectorGenome> children;
            std::size_t attempts=0;
            while(children.size()<evo.offspring_per_generation&&attempts++<50*evo.offspring_per_generation)
            {
                std::vector<std::size_t>idx(parents.size());
                std::iota(idx.begin(),idx.end(),0);
                std::shuffle(idx.begin(),idx.end(),rng);
                idx.resize(std::min<std::size_t>(3,idx.size()));
                auto best=idx.front();
                for(auto i:idx)if(key(parents[i])<key(parents[best]))best=i;
                auto child=mutate(parents[best].algorithm.selector,space,rng);
                if(std::find(children.begin(),children.end(),child)==children.end())children.push_back(child);
            }

            std::map<SelectorGenome,PortfolioResult>pool;
            for(auto&r:parents)pool[r.algorithm.selector]=r;
            for(auto&c:children)pool[c]=score(c);
            std::vector<PortfolioResult>ranked;
            for(auto&[_,r]:pool)ranked.push_back(r);
            std::sort(ranked.begin(),ranked.end(),[](auto&a,auto&b)
            {
                return key(a)<key(b);
            }

            );
            pop.clear();
            for(std::size_t i=0;i<std::min(evo.population_size,ranked.size());++i)pop.push_back(ranked[i].algorithm.selector);
            double best=std::numeric_limits<double>::infinity();
            for(auto&r:ranked)if(r.is_exact())best=std::min(best,r.mean_proxy_cost);
            history.push_back(best);
        }

        std::vector<PortfolioResult> evaluated;
        for(auto&[_,r]:cache)evaluated.push_back(r);
        std::sort(evaluated.begin(),evaluated.end(),[](auto&a,auto&b)
        {
            return key(a)<key(b);
        }

        );
        auto it=std::find_if(evaluated.begin(),evaluated.end(),[](auto&r)
        {
            return r.is_exact();
        }

        );
        if(it==evaluated.end())throw std::runtime_error("selector search found no exact algorithm");
        return
        {
            evaluated,*it,history
        };
    }
}

// namespace polar_scl
