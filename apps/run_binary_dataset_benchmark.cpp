#include "polar_scl/polar_scl.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <optional>
#include <string>
#include <vector>

using namespace polar_scl;

int main(int argc,char** argv)
{
    try
    {
        std::filesystem::path data;
        std::filesystem::path output="results/mass_validation/benchmark.csv";
        std::optional<std::filesystem::path> selector_path;
        std::size_t list_size=8;

        for(int index=1;index<argc;++index)
        {
            const std::string argument=argv[index];
            if(argument=="--data"&&index+1<argc)
            {
                data=argv[++index];
            }
            else if(argument=="--out"&&index+1<argc)
            {
                output=argv[++index];
            }
            else if(argument=="--selector"&&index+1<argc)
            {
                selector_path=std::filesystem::path(argv[++index]);
            }
            else if(argument=="--list"&&index+1<argc)
            {
                list_size=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else
            {
                throw std::invalid_argument("unknown argument: "+argument);
            }
        }
        if(data.empty())
        {
            throw std::invalid_argument("--data is required");
        }

        PolarAlgorithmIR searched_algorithm;
        if(selector_path)
        {
            searched_algorithm.selector_mode=SelectorMode::Rule;
            searched_algorithm.selector=load_selector_genome(*selector_path);
        }

        const auto dataset=load_dataset_binary(data);
        ExperimentRunConfig run_config;
        run_config.list_size=list_size;
        run_config.compare_with_scl=true;
        run_config.warmup_frames=std::min<std::size_t>(16,dataset.frames.size());
        run_config.timing_repeats=1;

        const auto variants=paper_decoder_variants(selector_path.has_value());
        std::vector<ExperimentMetrics> metrics;
        metrics.reserve(variants.size());
        for(DecoderVariant variant:variants)
        {
            const PolarAlgorithmIR algorithm=variant==DecoderVariant::SearchedPortfolio
                ?searched_algorithm:PolarAlgorithmIR{};
            const auto result=run_decoder_experiment(dataset,variant,run_config,algorithm);
            metrics.push_back(result);
            std::cout<<"[BENCH] "<<to_string(variant)
                <<" FER="<<result.fer()
                <<" mean_us="<<result.mean_us
                <<" p99_us="<<result.p99_us
                <<" mismatch="<<result.list_mismatches<<'\n';
        }
        write_metrics_csv(metrics,output);
        std::cout<<"benchmark CSV: "<<output.string()<<'\n';
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"binary benchmark failed: "<<error.what()<<'\n';
        return 1;
    }
}
