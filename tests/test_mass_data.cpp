#include "polar_scl/polar_scl.hpp"
#include "test_util.hpp"

#include <filesystem>

using namespace polar_scl;

int main()
{
    const auto temporary=std::filesystem::temp_directory_path()/"sca_polar_mass_data_test";
    std::filesystem::create_directories(temporary);

    AwgnDatasetConfig awgn;
    awgn.length=32;
    awgn.payload_bits=8;
    awgn.crc_bits=8;
    awgn.ebn0_db=1.0;
    awgn.frames=12;
    awgn.seed=12345;
    const auto dataset=generate_awgn_dataset(awgn);
    const auto dataset_path=temporary/"dataset.scabin";
    save_dataset_binary(dataset,dataset_path);
    const auto loaded=load_dataset_binary(dataset_path);
    CHECK_TRUE(loaded.config.length==dataset.config.length);
    CHECK_TRUE(loaded.frames.size()==dataset.frames.size());
    CHECK_TRUE(loaded.frames.front().payload==dataset.frames.front().payload);
    CHECK_TRUE(loaded.frames.front().llrs==dataset.frames.front().llrs);

    StressDatasetConfig stress_config;
    stress_config.cases_per_mode=4;
    stress_config.seed=67890;
    stress_config.lengths={32};
    stress_config.list_sizes={2,4};
    const auto stress=generate_stress_dataset(stress_config);
    CHECK_TRUE(stress.cases.size()==76);
    const auto stress_path=temporary/"stress.scastress";
    save_stress_dataset_binary(stress,stress_path);
    const auto loaded_stress=load_stress_dataset_binary(stress_path);
    CHECK_TRUE(loaded_stress.cases.size()==stress.cases.size());
    CHECK_TRUE(loaded_stress.cases.front().llrs==stress.cases.front().llrs);
    CHECK_TRUE(loaded_stress.cases.front().frozen==stress.cases.front().frozen);

    const auto report=run_stress_validation(loaded_stress,PolarAlgorithmIR{},false);
    std::size_t guarded_mismatches=0;
    std::size_t direct_mismatches=0;
    std::size_t raw_mismatches=0;
    for(const auto& row:report.rows)
    {
        guarded_mismatches+=row.guarded_mismatches;
        direct_mismatches+=row.direct_mismatches;
        raw_mismatches+=row.raw_mismatches;
    }
    CHECK_TRUE(guarded_mismatches==0);
    CHECK_TRUE(direct_mismatches==0);
    CHECK_TRUE(raw_mismatches>0);


    const std::vector<bool> all_information_128(128,false);
    const auto large_schedule=rate1_schedule(all_information_128);
    CHECK_EQ(large_schedule.size(),std::size_t{2});
    CHECK_TRUE(large_schedule[0].kind==DecodeOp::Kind::Rate1);
    CHECK_TRUE(large_schedule[1].kind==DecodeOp::Kind::Rate1);
    CHECK_EQ(large_schedule[0].size,std::size_t{64});
    CHECK_EQ(large_schedule[1].size,std::size_t{64});

    StressDatasetConfig boundary_config;
    boundary_config.cases_per_mode=2;
    boundary_config.seed=24680;
    boundary_config.lengths={128};
    boundary_config.list_sizes={8,16};
    const auto boundary_stress=generate_stress_dataset(boundary_config);
    const auto boundary_report=run_stress_validation(boundary_stress,PolarAlgorithmIR{},false);
    std::size_t boundary_mismatches=0;
    for(const auto& row:boundary_report.rows)
    {
        boundary_mismatches+=row.guarded_mismatches;
        boundary_mismatches+=row.direct_mismatches;
    }
    CHECK_EQ(boundary_mismatches,std::size_t{0});

    std::filesystem::remove_all(temporary);
    return 0;
}
