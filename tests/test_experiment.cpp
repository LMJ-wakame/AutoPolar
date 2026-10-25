#include "test_util.hpp"

#include <filesystem>

using namespace polar_scl;

int main()
{
    return run_test_main("experiment",[]
    {
        AwgnDatasetConfig config;
        config.length=32;
        config.payload_bits=8;
        config.crc_bits=8;
        config.ebn0_db=1.5;
        config.frames=12;
        config.seed=123456;

        const auto first=generate_awgn_dataset(config);
        const auto second=generate_awgn_dataset(config);
        CHECK_EQ(first.frames.size(),second.frames.size());
        for(std::size_t index=0;index<first.frames.size();++index)
        {
            CHECK_EQ(first.frames[index].payload,second.frames[index].payload);
            CHECK_EQ(first.frames[index].llrs,second.frames[index].llrs);
        }

        const auto path=std::filesystem::temp_directory_path()/"polar_scl_dataset_roundtrip.csv";
        save_dataset_csv(first,path);
        const auto loaded=load_dataset_csv(path);
        CHECK_EQ(loaded.config.length,config.length);
        CHECK_EQ(loaded.config.payload_bits,config.payload_bits);
        CHECK_EQ(loaded.config.crc_bits,config.crc_bits);
        CHECK_EQ(loaded.config.frames,config.frames);
        CHECK_EQ(loaded.frames.size(),first.frames.size());
        for(std::size_t index=0;index<first.frames.size();++index)
        {
            CHECK_EQ(loaded.frames[index].payload,first.frames[index].payload);
            CHECK_EQ(loaded.frames[index].llrs,first.frames[index].llrs);
        }

        std::filesystem::remove(path);

        ExperimentRunConfig run;
        run.list_size=4;
        run.compare_with_scl=true;
        run.warmup_frames=2;
        const auto metrics=run_decoder_experiment(first,DecoderVariant::CertifiedFrontier,run);
        CHECK_EQ(metrics.frames,config.frames);
        CHECK_EQ(metrics.compared_frames,config.frames);
        CHECK_EQ(metrics.list_mismatches,0U);
        CHECK_TRUE(metrics.mean_us>0.0);
        CHECK_TRUE(metrics.p99_us>=0.0);
    }
    );
}
