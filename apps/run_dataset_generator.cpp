#include "polar_scl/polar_scl.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace polar_scl;

namespace
{
    std::size_t parse_frames(int argc,char** argv)
    {
        std::size_t frames=128;
        for(int index=1;index<argc;++index)
        {
            const std::string argument=argv[index];
            if(argument=="--frames"&&index+1<argc)
            {
                frames=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
        }

        return frames;
    }

    std::filesystem::path parse_output(int argc,char** argv)
    {
        std::filesystem::path output="datasets/generated";
        for(int index=1;index<argc;++index)
        {
            const std::string argument=argv[index];
            if(argument=="--out"&&index+1<argc)
            {
                output=argv[++index];
            }
        }

        return output;
    }
}
// namespace

int main(int argc,char** argv)
{
    try
    {
        const std::size_t frames=parse_frames(argc,argv);
        const auto output_directory=parse_output(argc,argv);
        std::filesystem::create_directories(output_directory);

        std::vector<AwgnDatasetConfig> configs;
        for(double ebn0:{0.0,1.0,2.0})
        {
            configs.push_back({64,32,8,ebn0,frames,2026081600ULL+static_cast<std::uint64_t>(ebn0*10.0)});
        }

        for(double ebn0:{1.0,2.0})
        {
            configs.push_back({128,56,8,ebn0,frames,2026081700ULL+static_cast<std::uint64_t>(ebn0*10.0)});
        }

        std::cout<<"Generating reproducible AWGN datasets\n";
        std::cout<<"output="<<output_directory.string()<<" frames/config="<<frames<<'\n';
        for(const auto& config:configs)
        {
            const auto dataset=generate_awgn_dataset(config);
            std::ostringstream filename;
            filename<<"awgn_N"<<config.length
            <<"_K"<<config.information_bits()
            <<"_EbN0_"<<std::fixed<<std::setprecision(1)<<config.ebn0_db
            <<"dB.csv";
            const auto path=output_directory/filename.str();
            save_dataset_csv(dataset,path);
            std::cout<<"[DATASET] "<<path.string()
            <<" rate="<<std::setprecision(3)<<config.code_rate()
            <<" frames="<<dataset.frames.size()<<'\n';
        }

        return 0;
    }

    catch(const std::exception& error)
    {
        std::cerr<<"dataset generation failed: "<<error.what()<<'\n';
        return 1;
    }
}
