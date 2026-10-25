#include "polar_scl/polar_scl.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace polar_scl;

namespace
{
    struct Options
    {
        std::string profile{"quick"};
        std::size_t train_frames{12};
        std::size_t validation_frames{24};
        std::size_t population{10};
        std::size_t offspring{10};
        std::size_t generations{6};
        std::uint64_t seed{20260816};
        std::string objective{"latency"};
        std::size_t timing_repeats{2};
        std::filesystem::path output_directory{"results/sca_search"};
    };

    Options parse_options(int argc,char** argv)
    {
        Options options;
        for(int index=1;index<argc;++index)
        {
            const std::string argument=argv[index];
            if(argument=="--profile"&&index+1<argc)
            {
                options.profile=argv[++index];
            }
            else if(argument=="--train-frames"&&index+1<argc)
            {
                options.train_frames=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--validation-frames"&&index+1<argc)
            {
                options.validation_frames=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--population"&&index+1<argc)
            {
                options.population=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--offspring"&&index+1<argc)
            {
                options.offspring=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--generations"&&index+1<argc)
            {
                options.generations=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--seed"&&index+1<argc)
            {
                options.seed=static_cast<std::uint64_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--objective"&&index+1<argc)
            {
                options.objective=argv[++index];
            }
            else if(argument=="--timing-repeats"&&index+1<argc)
            {
                options.timing_repeats=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--out"&&index+1<argc)
            {
                options.output_directory=argv[++index];
            }
            else
            {
                throw std::invalid_argument("unknown argument: "+argument);
            }
        }

        if(options.profile!="quick"&&options.profile!="paper")
        {
            throw std::invalid_argument("profile must be quick or paper");
        }

        if(options.train_frames==0||options.validation_frames==0||options.timing_repeats==0)
        {
            throw std::invalid_argument("train/validation frames and timing repeats must be positive");
        }

        if(options.objective!="latency"&&options.objective!="proxy")
        {
            throw std::invalid_argument("objective must be latency or proxy");
        }

        return options;
    }

    std::vector<std::size_t> lengths(const Options& options)
    {
        return options.profile=="quick"?std::vector<std::size_t>{64,128}
            :std::vector<std::size_t>{128,256,512};
    }

    std::vector<double> rates()
    {
        return {0.5,0.75};
    }

    std::vector<double> ebn0_values(const Options& options)
    {
        return options.profile=="quick"?std::vector<double>{0.5,1.5,2.5}
            :std::vector<double>{0.5,1.5,2.5};
    }

    std::vector<std::size_t> list_sizes(const Options& options)
    {
        return options.profile=="quick"?std::vector<std::size_t>{4,8}
            :std::vector<std::size_t>{4,8,16};
    }

    std::uint64_t scenario_seed(
        std::uint64_t base,
        std::uint64_t split_offset,
        std::size_t length,
        std::size_t information_bits,
        double ebn0)
    {
        return base+split_offset
            +static_cast<std::uint64_t>(length)*100000ULL
            +static_cast<std::uint64_t>(information_bits)*100ULL
            +static_cast<std::uint64_t>(std::llround(ebn0*10.0));
    }

    DatasetSelectorEvaluationConfig make_evaluation(
        const Options& options,
        std::size_t frames,
        std::uint64_t split_offset)
    {
        DatasetSelectorEvaluationConfig evaluation;
        constexpr std::size_t crc_bits=8;
        for(std::size_t length:lengths(options))
        {
            for(double rate:rates())
            {
                const std::size_t information_bits=static_cast<std::size_t>(
                    std::llround(rate*static_cast<double>(length)));
                const std::size_t payload_bits=information_bits-crc_bits;
                for(double ebn0:ebn0_values(options))
                {
                    AwgnDatasetConfig config;
                    config.length=length;
                    config.payload_bits=payload_bits;
                    config.crc_bits=crc_bits;
                    config.ebn0_db=ebn0;
                    config.frames=frames;
                    config.seed=scenario_seed(
                        options.seed,split_offset,length,information_bits,ebn0);
                    evaluation.scenarios.push_back({generate_awgn_dataset(config),list_sizes(options)});
                }
            }
        }

        evaluation.objective=options.objective=="latency"
            ?SelectorObjective::Latency:SelectorObjective::ComplexityProxy;
        evaluation.timing_repeats=options.timing_repeats;
        return evaluation;
    }

    void write_summary(const SCASearchResult& result,const Options& options)
    {
        const auto path=options.output_directory/"SUMMARY.md";
        std::ofstream output(path);
        const auto& train=result.train_best;
        const auto& validation=result.validation_best;
        const auto& genome=validation.algorithm.selector;
        output<<"# SCA-Polar selector search summary\n\n";
        output<<"Search profile: `"<<options.profile<<"`  \n";
        output<<"Train frames/config: "<<options.train_frames<<"  \n";
        output<<"Validation frames/config: "<<options.validation_frames<<"  \n";
        output<<"Search objective: `"<<options.objective<<"`  \n";
        output<<"Timing repeats: "<<options.timing_repeats<<"  \n";
        output<<"Train and validation use disjoint random seeds. Final paper tests must use a third seed range.\n\n";
        output<<"## Best train candidate\n\n";
        output<<"- exact: "<<train.exact_cases<<"/"<<train.total_cases<<"\n";
        output<<"- mean proxy: "<<train.mean_proxy_cost<<"\n";
        output<<"- p99 proxy: "<<train.p99_proxy_cost<<"\n";
        output<<"- mean latency (search machine): "<<train.mean_latency_us<<" us\n";
        output<<"- objective cost: "<<train.objective_cost<<"\n\n";
        output<<"## Selected validation candidate (SCA-Polar)\n\n";
        output<<"- exact: "<<validation.exact_cases<<"/"<<validation.total_cases<<"\n";
        output<<"- mean proxy: "<<validation.mean_proxy_cost<<"\n";
        output<<"- p99 proxy: "<<validation.p99_proxy_cost<<"\n";
        output<<"- mean latency (search machine): "<<validation.mean_latency_us<<" us\n";
        output<<"- objective cost: "<<validation.objective_cost<<"\n";
        output<<"- list size <= "<<genome.list_size_threshold<<": "<<to_string(genome.low_list_engine)<<"\n";
        output<<"- small node <= "<<genome.small_node_threshold<<": "<<to_string(genome.small_engine)<<"\n";
        output<<"- easy engine: "<<to_string(genome.easy_engine)<<"\n";
        output<<"- hard engine: "<<to_string(genome.hard_engine)<<"\n";
        output<<"- reliability gap threshold: "<<genome.reliability_gap_threshold<<"\n";
        output<<"- parent metric gap threshold: "<<genome.parent_metric_gap_threshold<<"\n";
        output<<"- frontier initial bits: "<<genome.frontier_initial_active_bits<<"\n";
        output<<"- frontier growth step: "<<genome.frontier_growth_step<<"\n";
        output<<"- frontier budget: "<<genome.frontier_candidate_budget<<"\n";
        output<<"- REP special node: "<<(genome.enable_repetition?"enabled":"disabled")<<"\n";
        output<<"- SPC special node: "<<(genome.enable_spc?"enabled":"disabled")<<"\n\n";
        if(options.objective=="proxy")
        {
            output<<"The search objective is a deterministic decoding-work proxy.\n";
        }
        else
        {
            output<<"The search objective is measured C++ wall-clock latency on the search machine. Validation and final paper tests use disjoint frames.\n";
        }
    }
}

int main(int argc,char** argv)
{
    try
    {
        const Options options=parse_options(argc,argv);
        std::filesystem::create_directories(options.output_directory);
        std::cout<<"SCA-Polar search profile="<<options.profile
            <<" train_frames/config="<<options.train_frames
            <<" validation_frames/config="<<options.validation_frames
            <<" objective="<<options.objective<<'\n';

        const DatasetSelectorEvaluationConfig train=make_evaluation(
            options,options.train_frames,100000000ULL);
        const DatasetSelectorEvaluationConfig validation=make_evaluation(
            options,options.validation_frames,200000000ULL);

        SCASearchConfig search;
        search.evolution.population_size=options.population;
        search.evolution.offspring_per_generation=options.offspring;
        search.evolution.generations=options.generations;
        search.evolution.seed=options.seed;
        search.validation_top_k=std::max<std::size_t>(8,options.population);
        const SCASearchResult result=search_sca_selector(train,validation,search);

        write_selector_search_csv(
            result.train_evaluated,options.output_directory/"train_candidates.csv");
        write_selector_search_csv(
            result.validation_evaluated,options.output_directory/"validation_candidates.csv");
        save_selector_genome(
            result.validation_best.algorithm.selector,options.output_directory/"best_selector.txt");
        write_summary(result,options);

        std::cout<<"Best train objective="<<result.train_best.objective_cost
            <<" proxy="<<result.train_best.mean_proxy_cost
            <<" latency_us="<<result.train_best.mean_latency_us
            <<" exact="<<result.train_best.exact_cases<<'/'<<result.train_best.total_cases<<'\n';
        std::cout<<"Selected validation objective="<<result.validation_best.objective_cost
            <<" proxy="<<result.validation_best.mean_proxy_cost
            <<" latency_us="<<result.validation_best.mean_latency_us
            <<" exact="<<result.validation_best.exact_cases<<'/'<<result.validation_best.total_cases<<'\n';
        const SelectorGenome& genome=result.validation_best.algorithm.selector;
        std::cout<<"SCA-Polar selector: L<= "<<genome.list_size_threshold
            <<" -> "<<to_string(genome.low_list_engine)
            <<", small<= "<<genome.small_node_threshold
            <<" -> "<<to_string(genome.small_engine)
            <<", easy -> "<<to_string(genome.easy_engine)
            <<", hard -> "<<to_string(genome.hard_engine)<<'\n';
        std::cout<<"Saved selector: "<<(options.output_directory/"best_selector.txt").string()<<'\n';
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"[ERROR] "<<error.what()<<'\n';
        return 1;
    }
}
