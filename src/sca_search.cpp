#include "polar_scl/sca_search.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace polar_scl
{
    namespace
    {
        template<class T>
        T choose(const std::vector<T>& values,std::mt19937_64& rng)
        {
            if(values.empty())
            {
                throw std::invalid_argument("search-space slot must not be empty");
            }

            std::uniform_int_distribution<std::size_t> distribution(0,values.size()-1);
            return values[distribution(rng)];
        }

        template<class T>
        T choose_different(const std::vector<T>& values,const T& current,std::mt19937_64& rng)
        {
            std::vector<T> alternatives;
            for(const T& value:values)
            {
                if(value!=current)
                {
                    alternatives.push_back(value);
                }
            }

            return alternatives.empty()?current:choose(alternatives,rng);
        }

        SelectorGenome random_selector(const SelectorSearchSpace& space,std::mt19937_64& rng)
        {
            SelectorGenome genome;
            genome.list_size_threshold=choose(space.list_size_thresholds,rng);
            genome.low_list_engine=choose(space.low_list_engines,rng);
            genome.small_node_threshold=choose(space.small_node_thresholds,rng);
            genome.reliability_gap_threshold=choose(space.reliability_gap_thresholds,rng);
            genome.parent_metric_gap_threshold=choose(space.parent_metric_gap_thresholds,rng);
            genome.small_engine=choose(space.small_engines,rng);
            genome.easy_engine=choose(space.easy_engines,rng);
            genome.hard_engine=choose(space.hard_engines,rng);
            genome.frontier_initial_active_bits=choose(space.frontier_initial_bits,rng);
            genome.frontier_growth_step=choose(space.frontier_growth_steps,rng);
            genome.frontier_expansion_mode=choose(space.frontier_modes,rng);
            genome.frontier_candidate_budget=choose(space.frontier_budgets,rng);
            genome.enable_repetition=choose(space.repetition_choices,rng);
            genome.enable_spc=choose(space.spc_choices,rng);
            validate_selector_genome(genome);
            return genome;
        }

        SelectorGenome mutate_selector(
            const SelectorGenome& genome,
            const SelectorSearchSpace& space,
            std::mt19937_64& rng)
        {
            SelectorGenome child=genome;
            std::uniform_int_distribution<int> slot_distribution(0,13);
            switch(slot_distribution(rng))
            {
                case 0:
                child.list_size_threshold=choose_different(space.list_size_thresholds,genome.list_size_threshold,rng);
                break;
                case 1:
                child.low_list_engine=choose_different(space.low_list_engines,genome.low_list_engine,rng);
                break;
                case 2:
                child.small_node_threshold=choose_different(space.small_node_thresholds,genome.small_node_threshold,rng);
                break;
                case 3:
                child.reliability_gap_threshold=choose_different(
                    space.reliability_gap_thresholds,genome.reliability_gap_threshold,rng);
                break;
                case 4:
                child.parent_metric_gap_threshold=choose_different(
                    space.parent_metric_gap_thresholds,genome.parent_metric_gap_threshold,rng);
                break;
                case 5:
                child.small_engine=choose_different(space.small_engines,genome.small_engine,rng);
                break;
                case 6:
                child.easy_engine=choose_different(space.easy_engines,genome.easy_engine,rng);
                break;
                case 7:
                child.hard_engine=choose_different(space.hard_engines,genome.hard_engine,rng);
                break;
                case 8:
                child.frontier_initial_active_bits=choose_different(
                    space.frontier_initial_bits,genome.frontier_initial_active_bits,rng);
                break;
                case 9:
                child.frontier_growth_step=choose_different(
                    space.frontier_growth_steps,genome.frontier_growth_step,rng);
                break;
                case 10:
                child.frontier_expansion_mode=choose_different(
                    space.frontier_modes,genome.frontier_expansion_mode,rng);
                break;
                case 11:
                child.frontier_candidate_budget=choose_different(
                    space.frontier_budgets,genome.frontier_candidate_budget,rng);
                break;
                case 12:
                child.enable_repetition=choose_different(space.repetition_choices,genome.enable_repetition,rng);
                break;
                case 13:
                child.enable_spc=choose_different(space.spc_choices,genome.enable_spc,rng);
                break;
            }

            validate_selector_genome(child);
            return child;
        }

        double percentile(std::vector<double> values,double probability)
        {
            if(values.empty())
            {
                return 0.0;
            }

            std::sort(values.begin(),values.end());
            const std::size_t index=static_cast<std::size_t>(
                std::floor(probability*static_cast<double>(values.size()-1)));
            return values[index];
        }

        double frame_proxy_cost(
            const PortfolioDecoderStats& stats,
            const FrameCostWeights& weights)
        {
            const DecoderStats& decoder=stats.decoder;
            return weights.leaf_step*static_cast<double>(decoder.leaf_steps)
                +weights.path_split*static_cast<double>(decoder.path_splits)
                +weights.prune_event*static_cast<double>(decoder.sorts)
                +weights.rate1_candidate*static_cast<double>(decoder.rate1_candidates)
                +weights.frontier_expansion*static_cast<double>(decoder.rate1_frontier_expansions)
                +weights.rate1_fallback*static_cast<double>(decoder.rate1_fallbacks)
                +weights.tie_fallback*static_cast<double>(stats.tie_fallbacks);
        }

        auto result_key(const DatasetSelectorResult& result)
        {
            return std::tuple{
                !result.is_exact(),
                result.objective_cost,
                result.objective_p99,
                result.fallback_rate,
                result.mean_candidates};
        }

        SelectorGenome fixed_selector(Rate1Engine engine)
        {
            SelectorGenome genome;
            genome.list_size_threshold=0;
            genome.small_node_threshold=2;
            genome.small_engine=engine;
            genome.easy_engine=engine;
            genome.hard_engine=engine;
            genome.frontier_initial_active_bits=0;
            genome.frontier_growth_step=1;
            genome.frontier_expansion_mode=ExpansionMode::MinimumBound;
            genome.frontier_candidate_budget=256;
            return genome;
        }

        std::vector<SelectorGenome> seed_genomes()
        {
            std::vector<SelectorGenome> seeds;
            seeds.reserve(4);
            seeds.push_back(fixed_selector(Rate1Engine::LazyExact));
            seeds.push_back(fixed_selector(Rate1Engine::FastSSCLBounded));
            seeds.push_back(fixed_selector(Rate1Engine::CertifiedFrontier));

            SelectorGenome rule;
            rule.list_size_threshold=4;
            rule.low_list_engine=Rate1Engine::CertifiedFrontier;
            rule.small_node_threshold=4;
            rule.small_engine=Rate1Engine::LazyExact;
            rule.easy_engine=Rate1Engine::CertifiedFrontier;
            rule.hard_engine=Rate1Engine::FastSSCLBounded;
            rule.reliability_gap_threshold=1.0;
            rule.parent_metric_gap_threshold=0.5;
            rule.frontier_initial_active_bits=0;
            rule.frontier_growth_step=1;
            rule.frontier_expansion_mode=ExpansionMode::MinimumBound;
            rule.frontier_candidate_budget=256;
            seeds.push_back(rule);
            return seeds;
        }

        std::string engine_token(Rate1Engine engine)
        {
            return to_string(engine);
        }

        Rate1Engine parse_engine(const std::string& text)
        {
            if(text=="exhaustive")
            {
                return Rate1Engine::Exhaustive;
            }

            if(text=="lazy_exact")
            {
                return Rate1Engine::LazyExact;
            }

            if(text=="fast_sscl_bounded")
            {
                return Rate1Engine::FastSSCLBounded;
            }

            if(text=="certified_frontier")
            {
                return Rate1Engine::CertifiedFrontier;
            }

            throw std::invalid_argument("unknown Rate-1 engine token: "+text);
        }

        std::string expansion_token(ExpansionMode mode)
        {
            return mode==ExpansionMode::MinimumBound?"minimum_bound":"all_competitive";
        }

        ExpansionMode parse_expansion(const std::string& text)
        {
            if(text=="minimum_bound")
            {
                return ExpansionMode::MinimumBound;
            }

            if(text=="all_competitive")
            {
                return ExpansionMode::AllCompetitive;
            }

            throw std::invalid_argument("unknown expansion mode token: "+text);
        }

        std::map<std::string,std::string> read_key_values(const std::filesystem::path& path)
        {
            std::ifstream input(path);
            if(!input)
            {
                throw std::runtime_error("failed to open selector file: "+path.string());
            }

            std::map<std::string,std::string> values;
            std::string line;
            while(std::getline(input,line))
            {
                if(line.empty()||line[0]=='#')
                {
                    continue;
                }

                const std::size_t position=line.find('=');
                if(position==std::string::npos)
                {
                    throw std::runtime_error("invalid selector line: "+line);
                }

                values[line.substr(0,position)]=line.substr(position+1);
            }

            return values;
        }
    }

    DatasetSelectorResult evaluate_selector_on_datasets(
        const SelectorGenome& genome,
        const DatasetSelectorEvaluationConfig& config)
    {
        validate_selector_genome(genome);
        if(config.scenarios.empty())
        {
            throw std::invalid_argument("at least one selector-search scenario is required");
        }

        if(config.timing_repeats==0)
        {
            throw std::invalid_argument("timing_repeats must be positive");
        }

        PolarAlgorithmIR algorithm;
        algorithm.selector_mode=SelectorMode::Rule;
        algorithm.selector=genome;

        std::vector<double> proxy_costs;
        std::vector<double> latencies;
        double candidate_sum=0.0;
        std::size_t exact_cases=0;
        std::size_t total_cases=0;
        std::size_t fallback_cases=0;
        std::map<Rate1Engine,std::size_t> engine_counts;

        for(const SelectorScenario& scenario:config.scenarios)
        {
            const Dataset& dataset=scenario.dataset;
            const auto frozen=pw_frozen_pattern(dataset.config.length,dataset.config.information_bits());
            if(scenario.list_sizes.empty())
            {
                throw std::invalid_argument("selector-search scenario has no list sizes");
            }

            for(std::size_t list_size:scenario.list_sizes)
            {
                if(config.objective==SelectorObjective::Latency&&!dataset.frames.empty())
                {
                    const auto warmup=portfolio_scl_decode(
                        dataset.frames.front().llrs,frozen,list_size,algorithm);
                    if(warmup.paths.empty())
                    {
                        throw std::runtime_error("selector-search warmup returned no paths");
                    }
                }

                for(const ExperimentFrame& frame:dataset.frames)
                {
                    const DecodeResult reference=scl_decode(frame.llrs,frozen,list_size);
                    PortfolioDecodeResult actual;
                    for(std::size_t repeat=0;repeat<config.timing_repeats;++repeat)
                    {
                        const auto start_time=std::chrono::steady_clock::now();
                        PortfolioDecodeResult timed=portfolio_scl_decode(
                            frame.llrs,frozen,list_size,algorithm);
                        const auto stop_time=std::chrono::steady_clock::now();
                        if(config.objective==SelectorObjective::Latency)
                        {
                            latencies.push_back(
                                std::chrono::duration<double,std::micro>(stop_time-start_time).count());
                        }

                        if(repeat==0)
                        {
                            actual=std::move(timed);
                        }
                    }

                    exact_cases+=paths_equivalent(reference.paths,actual.paths)?1U:0U;
                    ++total_cases;
                    proxy_costs.push_back(frame_proxy_cost(actual.stats,config.weights));
                    candidate_sum+=static_cast<double>(actual.stats.decoder.rate1_candidates);
                    fallback_cases+=(actual.stats.decoder.rate1_fallbacks>0||actual.stats.tie_fallbacks>0)?1U:0U;
                    for(const auto& [engine,count]:actual.stats.engine_counts)
                    {
                        engine_counts[engine]+=count;
                    }
                }
            }
        }

        const double mean_proxy=std::accumulate(proxy_costs.begin(),proxy_costs.end(),0.0)
            /static_cast<double>(proxy_costs.size());
        const double mean_latency=latencies.empty()?0.0:
            std::accumulate(latencies.begin(),latencies.end(),0.0)/static_cast<double>(latencies.size());
        const double p99_proxy=percentile(proxy_costs,0.99);
        const double p99_latency=percentile(latencies,0.99);
        const bool latency_objective=config.objective==SelectorObjective::Latency;
        std::vector<std::pair<Rate1Engine,std::size_t>> counts(engine_counts.begin(),engine_counts.end());

        DatasetSelectorResult result;
        result.algorithm=algorithm;
        result.exact_cases=exact_cases;
        result.total_cases=total_cases;
        result.mean_proxy_cost=mean_proxy;
        result.p99_proxy_cost=p99_proxy;
        result.mean_candidates=candidate_sum/static_cast<double>(total_cases);
        result.fallback_rate=static_cast<double>(fallback_cases)/static_cast<double>(total_cases);
        result.mean_latency_us=mean_latency;
        result.p99_latency_us=p99_latency;
        result.objective_cost=latency_objective?mean_latency:mean_proxy;
        result.objective_p99=latency_objective?p99_latency:p99_proxy;
        result.engine_counts=std::move(counts);
        return result;
    }

    SCASearchResult search_sca_selector(
        const DatasetSelectorEvaluationConfig& train,
        const DatasetSelectorEvaluationConfig& validation,
        const SCASearchConfig& config)
    {
        if(config.evolution.population_size==0||
            config.evolution.offspring_per_generation==0||
            config.evolution.generations==0)
        {
            throw std::invalid_argument("population, offspring and generations must be positive");
        }

        if(config.validation_top_k==0)
        {
            throw std::invalid_argument("validation_top_k must be positive");
        }

        std::mt19937_64 rng(config.evolution.seed);
        std::map<SelectorGenome,DatasetSelectorResult> train_cache;
        auto score_train=[&](const SelectorGenome& genome)->const DatasetSelectorResult&
        {
            auto [iterator,inserted]=train_cache.try_emplace(genome);
            if(inserted)
            {
                iterator->second=evaluate_selector_on_datasets(genome,train);
            }

            return iterator->second;
        };

        std::vector<SelectorGenome> population;
        for(const SelectorGenome& seed:seed_genomes())
        {
            if(population.size()>=config.evolution.population_size)
            {
                break;
            }

            if(std::find(population.begin(),population.end(),seed)==population.end())
            {
                population.push_back(seed);
            }
        }

        while(population.size()<config.evolution.population_size)
        {
            SelectorGenome genome=random_selector(config.search_space,rng);
            if(std::find(population.begin(),population.end(),genome)==population.end())
            {
                population.push_back(genome);
            }
        }

        std::vector<double> history;
        for(std::size_t generation=0;generation<config.evolution.generations;++generation)
        {
            std::vector<DatasetSelectorResult> parents;
            for(const SelectorGenome& genome:population)
            {
                parents.push_back(score_train(genome));
            }

            std::sort(parents.begin(),parents.end(),[](const auto& lhs,const auto& rhs)
            {
                return result_key(lhs)<result_key(rhs);
            });

            std::vector<SelectorGenome> children;
            std::size_t attempts=0;
            while(children.size()<config.evolution.offspring_per_generation&&
                attempts<100*config.evolution.offspring_per_generation)
            {
                ++attempts;
                std::uniform_int_distribution<std::size_t> parent_distribution(
                    0,std::min<std::size_t>(parents.size(),4)-1);
                const SelectorGenome& parent=parents[parent_distribution(rng)].algorithm.selector;
                SelectorGenome child=mutate_selector(parent,config.search_space,rng);
                if(std::find(children.begin(),children.end(),child)==children.end())
                {
                    children.push_back(child);
                }
            }

            std::map<SelectorGenome,DatasetSelectorResult> pool;
            for(const DatasetSelectorResult& result:parents)
            {
                pool[result.algorithm.selector]=result;
            }

            for(const SelectorGenome& child:children)
            {
                pool[child]=score_train(child);
            }

            std::vector<DatasetSelectorResult> ranked;
            for(const auto& [genome,result]:pool)
            {
                (void)genome;
                ranked.push_back(result);
            }

            std::sort(ranked.begin(),ranked.end(),[](const auto& lhs,const auto& rhs)
            {
                return result_key(lhs)<result_key(rhs);
            });

            population.clear();
            for(std::size_t index=0;
                index<std::min(config.evolution.population_size,ranked.size());
                ++index)
            {
                population.push_back(ranked[index].algorithm.selector);
            }

            const auto best=std::find_if(ranked.begin(),ranked.end(),[](const auto& result)
            {
                return result.is_exact();
            });
            history.push_back(best==ranked.end()?std::numeric_limits<double>::infinity():best->objective_cost);
        }

        std::vector<DatasetSelectorResult> train_results;
        for(const auto& [genome,result]:train_cache)
        {
            (void)genome;
            train_results.push_back(result);
        }

        std::sort(train_results.begin(),train_results.end(),[](const auto& lhs,const auto& rhs)
        {
            return result_key(lhs)<result_key(rhs);
        });

        const auto train_best_iterator=std::find_if(
            train_results.begin(),train_results.end(),[](const auto& result)
        {
            return result.is_exact();
        });
        if(train_best_iterator==train_results.end())
        {
            throw std::runtime_error("SCA search found no exact training candidate");
        }

        std::vector<DatasetSelectorResult> validation_results;
        for(const DatasetSelectorResult& train_result:train_results)
        {
            if(!train_result.is_exact())
            {
                continue;
            }

            validation_results.push_back(
                evaluate_selector_on_datasets(train_result.algorithm.selector,validation));
            if(validation_results.size()>=config.validation_top_k)
            {
                break;
            }
        }

        std::sort(validation_results.begin(),validation_results.end(),[](const auto& lhs,const auto& rhs)
        {
            return result_key(lhs)<result_key(rhs);
        });

        const auto validation_best_iterator=std::find_if(
            validation_results.begin(),validation_results.end(),[](const auto& result)
        {
            return result.is_exact();
        });
        if(validation_best_iterator==validation_results.end())
        {
            throw std::runtime_error("SCA search found no exact validation candidate");
        }

        return {
            std::move(train_results),
            std::move(validation_results),
            *train_best_iterator,
            *validation_best_iterator,
            std::move(history)};
    }

    void save_selector_genome(
        const SelectorGenome& genome,
        const std::filesystem::path& path)
    {
        validate_selector_genome(genome);
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }

        std::ofstream output(path);
        if(!output)
        {
            throw std::runtime_error("failed to create selector file: "+path.string());
        }

        output<<std::setprecision(17);
        output<<"list_size_threshold="<<genome.list_size_threshold<<'\n';
        output<<"low_list_engine="<<engine_token(genome.low_list_engine)<<'\n';
        output<<"small_node_threshold="<<genome.small_node_threshold<<'\n';
        output<<"reliability_gap_threshold="<<genome.reliability_gap_threshold<<'\n';
        output<<"parent_metric_gap_threshold="<<genome.parent_metric_gap_threshold<<'\n';
        output<<"small_engine="<<engine_token(genome.small_engine)<<'\n';
        output<<"easy_engine="<<engine_token(genome.easy_engine)<<'\n';
        output<<"hard_engine="<<engine_token(genome.hard_engine)<<'\n';
        output<<"frontier_initial_active_bits="<<genome.frontier_initial_active_bits<<'\n';
        output<<"frontier_growth_step="<<genome.frontier_growth_step<<'\n';
        output<<"frontier_expansion_mode="<<expansion_token(genome.frontier_expansion_mode)<<'\n';
        output<<"frontier_candidate_budget="<<genome.frontier_candidate_budget<<'\n';
        output<<"enable_repetition="<<(genome.enable_repetition?1:0)<<'\n';
        output<<"enable_spc="<<(genome.enable_spc?1:0)<<'\n';
    }

    SelectorGenome load_selector_genome(const std::filesystem::path& path)
    {
        const auto values=read_key_values(path);
        const auto require=[&](const std::string& key)->const std::string&
        {
            const auto iterator=values.find(key);
            if(iterator==values.end())
            {
                throw std::runtime_error("selector file misses key: "+key);
            }

            return iterator->second;
        };

        SelectorGenome genome;
        genome.list_size_threshold=static_cast<std::size_t>(std::stoull(require("list_size_threshold")));
        genome.low_list_engine=parse_engine(require("low_list_engine"));
        genome.small_node_threshold=static_cast<std::size_t>(std::stoull(require("small_node_threshold")));
        genome.reliability_gap_threshold=std::stod(require("reliability_gap_threshold"));
        genome.parent_metric_gap_threshold=std::stod(require("parent_metric_gap_threshold"));
        genome.small_engine=parse_engine(require("small_engine"));
        genome.easy_engine=parse_engine(require("easy_engine"));
        genome.hard_engine=parse_engine(require("hard_engine"));
        genome.frontier_initial_active_bits=static_cast<std::size_t>(
            std::stoull(require("frontier_initial_active_bits")));
        genome.frontier_growth_step=static_cast<std::size_t>(std::stoull(require("frontier_growth_step")));
        genome.frontier_expansion_mode=parse_expansion(require("frontier_expansion_mode"));
        genome.frontier_candidate_budget=static_cast<std::size_t>(
            std::stoull(require("frontier_candidate_budget")));
        const auto repetition_iterator=values.find("enable_repetition");
        if(repetition_iterator!=values.end())
        {
            genome.enable_repetition=std::stoi(repetition_iterator->second)!=0;
        }
        const auto spc_iterator=values.find("enable_spc");
        if(spc_iterator!=values.end())
        {
            genome.enable_spc=std::stoi(spc_iterator->second)!=0;
        }
        validate_selector_genome(genome);
        return genome;
    }

    void write_selector_search_csv(
        const std::vector<DatasetSelectorResult>& results,
        const std::filesystem::path& path)
    {
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }

        std::ofstream output(path);
        if(!output)
        {
            throw std::runtime_error("failed to create selector-search CSV: "+path.string());
        }

        output<<"rank,exact_cases,total_cases,mean_proxy,p99_proxy,mean_latency_us,p99_latency_us,objective_cost,objective_p99,mean_candidates,fallback_rate,"
            "list_size_threshold,low_list_engine,small_threshold,reliability_gap_threshold,parent_metric_gap_threshold,small_engine,easy_engine,hard_engine,"
            "frontier_initial_bits,frontier_growth_step,frontier_mode,frontier_budget,enable_repetition,enable_spc\n";
        output<<std::setprecision(12);
        for(std::size_t index=0;index<results.size();++index)
        {
            const DatasetSelectorResult& result=results[index];
            const SelectorGenome& genome=result.algorithm.selector;
            output<<index+1<<','
                <<result.exact_cases<<','
                <<result.total_cases<<','
                <<result.mean_proxy_cost<<','
                <<result.p99_proxy_cost<<','
                <<result.mean_latency_us<<','
                <<result.p99_latency_us<<','
                <<result.objective_cost<<','
                <<result.objective_p99<<','
                <<result.mean_candidates<<','
                <<result.fallback_rate<<','
                <<genome.list_size_threshold<<','
                <<engine_token(genome.low_list_engine)<<','
                <<genome.small_node_threshold<<','
                <<genome.reliability_gap_threshold<<','
                <<genome.parent_metric_gap_threshold<<','
                <<engine_token(genome.small_engine)<<','
                <<engine_token(genome.easy_engine)<<','
                <<engine_token(genome.hard_engine)<<','
                <<genome.frontier_initial_active_bits<<','
                <<genome.frontier_growth_step<<','
                <<expansion_token(genome.frontier_expansion_mode)<<','
                <<genome.frontier_candidate_budget<<','
                <<(genome.enable_repetition?1:0)<<','
                <<(genome.enable_spc?1:0)<<'\n';
        }
    }
}
