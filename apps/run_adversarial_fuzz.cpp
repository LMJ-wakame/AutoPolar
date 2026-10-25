#include "polar_scl/polar_scl.hpp"

#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

using namespace polar_scl;

int main(int argc,char** argv)
{
    try
    {
        std::size_t cases_per_mode=1000;
        std::uint64_t seed=202608162222ULL;
        std::filesystem::path output="results/adversarial_fuzz";
        std::optional<std::filesystem::path> selector_path;

        for(int index=1;index<argc;++index)
        {
            const std::string argument=argv[index];
            if(argument=="--cases-per-mode"&&index+1<argc)
            {
                cases_per_mode=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--seed"&&index+1<argc)
            {
                seed=static_cast<std::uint64_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--out"&&index+1<argc)
            {
                output=argv[++index];
            }
            else if(argument=="--selector"&&index+1<argc)
            {
                selector_path=std::filesystem::path(argv[++index]);
            }
            else
            {
                throw std::invalid_argument("unknown argument: "+argument);
            }
        }

        StressDatasetConfig config;
        config.cases_per_mode=cases_per_mode;
        config.seed=seed;
        config.lengths={32,64,128};
        config.list_sizes={2,4,8,16};
        const auto dataset=generate_stress_dataset(config);

        PolarAlgorithmIR searched;
        if(selector_path)
        {
            searched.selector_mode=SelectorMode::Rule;
            searched.selector=load_selector_genome(*selector_path);
        }

        const auto report=run_stress_validation(dataset,searched,selector_path.has_value(),4096);
        std::filesystem::create_directories(output);
        const auto csv=output/"adversarial_validation.csv";
        const auto examples=output/"mismatch_examples.csv";
        const auto corpus=output/("adversarial_"+std::to_string(dataset.cases.size())+".scastress");
        write_stress_validation_csv(report.rows,csv);
        write_stress_mismatch_examples_csv(report.mismatch_examples,examples);
        save_stress_dataset_binary(dataset,corpus);

        std::size_t guarded=0;
        std::size_t direct=0;
        std::size_t raw=0;
        std::size_t fast_cases=0;
        for(const auto& row:report.rows)
        {
            guarded+=row.guarded_mismatches;
            direct+=row.direct_mismatches;
            raw+=row.raw_mismatches;
            fast_cases+=row.production_fast_path_cases();
        }

        std::cout<<"cases="<<dataset.cases.size()
            <<" guarded_mismatches="<<guarded
            <<" direct_mismatches="<<direct
            <<" raw_fastpath_mismatches="<<raw
            <<" production_fast_path_cases="<<fast_cases<<'\n';
        std::cout<<"CSV: "<<csv.string()<<'\n';
        std::cout<<"Examples: "<<examples.string()<<'\n';
        std::cout<<"Corpus: "<<corpus.string()<<'\n';
        return guarded==0&&direct==0?0:2;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"adversarial fuzz failed: "<<error.what()<<'\n';
        return 1;
    }
}
