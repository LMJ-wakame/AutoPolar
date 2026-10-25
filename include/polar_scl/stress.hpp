#pragma once

#include "polar_scl/algorithm_ir.hpp"
#include "polar_scl/decoder.hpp"
#include "polar_scl/experiment.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace polar_scl
{
    enum class StressMode:std::uint32_t
    {
        Continuous=0,
        Quantized8=1,
        Quantized6=2,
        Quantized4=3,
        Quantized3=4,
        NearZero=5,
        TieHeavy=6,
        MetricBoundary=7,
        RepetitionWholeNode=8,
        SpcWholeNode=9,
        RandomFrozen=10,
        SubsetSumCollision=11,
        GuardBoundaryOutside=12,
        ParentCompetition=13,
        DynamicRange=14,
        SignedZero=15,
        CancellationChain=16,
        DeepRandomFrozen=17,
        PermutedReliability=18
    };

    std::string to_string(StressMode mode);

    struct StressCase
    {
        StressMode mode{StressMode::Continuous};
        std::size_t list_size{8};
        std::vector<bool> frozen;
        std::vector<double> llrs;
    };

    struct StressDataset
    {
        std::uint64_t seed{20260816};
        std::vector<StressCase> cases;
    };

    struct StressDatasetConfig
    {
        std::size_t cases_per_mode{5000};
        std::uint64_t seed{20260816};
        std::vector<std::size_t> lengths{32,64,128};
        std::vector<std::size_t> list_sizes{2,4,8,16};
    };

    StressDataset generate_stress_dataset(const StressDatasetConfig& config);
    void save_stress_dataset_binary(const StressDataset& dataset,const std::filesystem::path& path);
    StressDataset load_stress_dataset_binary(const std::filesystem::path& path);

    struct StressValidationRow
    {
        DecoderVariant variant{DecoderVariant::SCL};
        StressMode mode{StressMode::Continuous};
        std::size_t cases{0};
        std::size_t outer_guard_cases{0};
        std::size_t guarded_mismatches{0};
        std::size_t guarded_tie_fallbacks{0};
        std::size_t guarded_rate1_fallbacks{0};
        std::size_t direct_mismatches{0};
        std::size_t direct_tie_fallbacks{0};
        std::size_t direct_rate1_fallbacks{0};
        std::size_t raw_mismatches{0};
        double guarded_mean_us{0.0};
        double direct_mean_us{0.0};
        double raw_mean_us{0.0};

        double guarded_mismatch_rate() const;
        double direct_mismatch_rate() const;
        double raw_mismatch_rate() const;
        std::size_t production_fast_path_cases() const;
    };

    enum class StressLayer
    {
        GuardedProduction,
        DirectInnerGuard,
        RawFastPath
    };

    std::string to_string(StressLayer layer);

    struct StressMismatchExample
    {
        DecoderVariant variant{DecoderVariant::SCL};
        StressMode mode{StressMode::Continuous};
        StressLayer layer{StressLayer::RawFastPath};
        std::size_t list_size{0};
        std::vector<bool> frozen;
        std::vector<double> llrs;
    };

    struct StressValidationReport
    {
        std::vector<StressValidationRow> rows;
        std::vector<StressMismatchExample> mismatch_examples;
    };

    StressValidationReport run_stress_validation(
        const StressDataset& dataset,
        const PolarAlgorithmIR& searched_algorithm={},
        bool include_searched=false,
        std::size_t max_saved_mismatches=256);

    void write_stress_validation_csv(
        const std::vector<StressValidationRow>& rows,
        const std::filesystem::path& path);

    void write_stress_mismatch_examples_csv(
        const std::vector<StressMismatchExample>& examples,
        const std::filesystem::path& path);
}
