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
        std::filesystem::path data="datasets/mass/stress_bin/adversarial_stress_95000.scastress";
        std::filesystem::path output="results/adversarial_validation/stress_validation.csv";
        std::filesystem::path examples="results/adversarial_validation/raw_mismatch_examples.csv";
        std::optional<std::filesystem::path> selector_path;

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
            else if(argument=="--examples"&&index+1<argc)
            {
                examples=argv[++index];
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

        PolarAlgorithmIR searched_algorithm;
        if(selector_path)
        {
            searched_algorithm.selector_mode=SelectorMode::Rule;
            searched_algorithm.selector=load_selector_genome(*selector_path);
        }

        const auto dataset=load_stress_dataset_binary(data);
        const auto report=run_stress_validation(
            dataset,
            searched_algorithm,
            selector_path.has_value(),
            5000);
        write_stress_validation_csv(report.rows,output);
        write_stress_mismatch_examples_csv(report.mismatch_examples,examples);

        std::size_t guarded_mismatches=0;
        std::size_t direct_mismatches=0;
        std::size_t raw_mismatches=0;
        std::size_t comparisons=0;
        std::size_t production_fast_path_cases=0;
        for(const auto& row:report.rows)
        {
            guarded_mismatches+=row.guarded_mismatches;
            direct_mismatches+=row.direct_mismatches;
            raw_mismatches+=row.raw_mismatches;
            comparisons+=row.cases;
            production_fast_path_cases+=row.production_fast_path_cases();
            std::cout<<"[STRESS] "<<to_string(row.variant)
                <<" mode="<<to_string(row.mode)
                <<" guarded="<<row.guarded_mismatches<<'/'<<row.cases
                <<" direct="<<row.direct_mismatches<<'/'<<row.cases
                <<" raw="<<row.raw_mismatches<<'/'<<row.cases
                <<" outer_guard="<<row.outer_guard_cases
                <<" fast_path="<<row.production_fast_path_cases()<<'\n';
        }
        std::cout<<"stress CSV: "<<output.string()<<'\n';
        std::cout<<"raw mismatch examples: "<<examples.string()<<'\n';
        std::cout<<"total comparisons/layer="<<comparisons
            <<" production_fast_path_cases="<<production_fast_path_cases
            <<" guarded_mismatches="<<guarded_mismatches
            <<" direct_mismatches="<<direct_mismatches
            <<" raw_fastpath_mismatches="<<raw_mismatches<<'\n';

        return guarded_mismatches==0&&direct_mismatches==0?0:2;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"stress validation failed: "<<error.what()<<'\n';
        return 1;
    }
}
