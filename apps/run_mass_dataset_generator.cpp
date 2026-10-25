#include "polar_scl/polar_scl.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace polar_scl;

namespace
{
    struct Options
    {
        std::filesystem::path output{"datasets/mass"};
        std::size_t stress_cases_per_mode{5000};
        bool stress_only{false};
    };

    Options parse_options(int argc,char** argv)
    {
        Options options;
        for(int index=1;index<argc;++index)
        {
            const std::string argument=argv[index];
            if(argument=="--out"&&index+1<argc)
            {
                options.output=argv[++index];
            }
            else if(argument=="--stress-cases-per-mode"&&index+1<argc)
            {
                options.stress_cases_per_mode=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--stress-only")
            {
                options.stress_only=true;
            }
            else
            {
                throw std::invalid_argument("unknown argument: "+argument);
            }
        }
        if(options.stress_cases_per_mode==0)
        {
            throw std::invalid_argument("stress-cases-per-mode must be positive");
        }
        return options;
    }

    std::string benchmark_name(const AwgnDatasetConfig& config)
    {
        std::ostringstream stream;
        stream<<"awgn_N"<<config.length
            <<"_K"<<config.information_bits()
            <<"_EbN0_"<<std::fixed<<std::setprecision(1)<<config.ebn0_db
            <<"_F"<<config.frames<<".scabin";
        return stream.str();
    }
}

int main(int argc,char** argv)
{
    try
    {
        const Options options=parse_options(argc,argv);
        const auto benchmark_directory=options.output/"benchmark_bin";
        const auto stress_directory=options.output/"stress_bin";
        std::filesystem::create_directories(stress_directory);

        std::size_t benchmark_frames=0;
        if(!options.stress_only)
        {
            std::filesystem::create_directories(benchmark_directory);
            std::vector<AwgnDatasetConfig> configs;
            std::uint64_t seed=202608160000ULL;
            for(double ebn0:{0.5,1.5,2.5})
            {
                configs.push_back({128,56,8,ebn0,1000,seed++});
                configs.push_back({128,88,8,ebn0,1000,seed++});
            }
            for(double ebn0:{1.0,2.0})
            {
                configs.push_back({256,120,8,ebn0,500,seed++});
                configs.push_back({256,184,8,ebn0,500,seed++});
            }
            configs.push_back({512,248,8,1.5,250,seed++});
            configs.push_back({512,376,8,1.5,250,seed++});

            std::cout<<"Generating mass benchmark datasets\n";
            for(const auto& config:configs)
            {
                const auto dataset=generate_awgn_dataset(config);
                const auto path=benchmark_directory/benchmark_name(config);
                save_dataset_binary(dataset,path);
                benchmark_frames+=dataset.frames.size();
                std::cout<<"[BENCH] "<<path.string()
                    <<" frames="<<dataset.frames.size()
                    <<" N="<<config.length
                    <<" K="<<config.information_bits()
                    <<" EbN0="<<config.ebn0_db<<'\n';
            }
        }

        StressDatasetConfig stress_config;
        stress_config.cases_per_mode=options.stress_cases_per_mode;
        stress_config.seed=202608169999ULL;
        stress_config.lengths={32,64,128};
        stress_config.list_sizes={2,4,8,16};
        const auto stress=generate_stress_dataset(stress_config);
        const auto stress_path=stress_directory/
            ("adversarial_stress_"+std::to_string(stress.cases.size())+".scastress");
        save_stress_dataset_binary(stress,stress_path);

        std::cout<<"[STRESS] "<<stress_path.string()
            <<" cases="<<stress.cases.size()<<'\n';
        std::cout<<"benchmark_frames="<<benchmark_frames
            <<" stress_cases="<<stress.cases.size()<<'\n';
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"mass dataset generation failed: "<<error.what()<<'\n';
        return 1;
    }
}
