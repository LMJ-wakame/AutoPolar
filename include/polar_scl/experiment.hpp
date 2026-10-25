#pragma once

#include "polar_scl/algorithm_ir.hpp"
#include "polar_scl/construction.hpp"
#include "polar_scl/decoder.hpp"
#include "polar_scl/portfolio.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace polar_scl
{
    enum class DecoderVariant
    {
        SC,
        SCL,
        SSCLRate01Exact,
        FastSSCLRate1Bounded,
        CertifiedFrontier,
        PortfolioRule,
        SearchedPortfolio
    };

    std::string to_string(DecoderVariant variant);
    std::vector<DecoderVariant> paper_decoder_variants(bool include_searched=false);

    struct AwgnDatasetConfig
    {
        std::size_t length{64};
        std::size_t payload_bits{32};
        std::size_t crc_bits{8};
        double ebn0_db{1.0};
        std::size_t frames{100};
        std::uint64_t seed{20260816};

        std::size_t information_bits() const;
        double code_rate() const;
    };

    struct ExperimentFrame
    {
        std::vector<Bit> payload;
        std::vector<double> llrs;
    };

    struct Dataset
    {
        AwgnDatasetConfig config;
        std::vector<ExperimentFrame> frames;
    };

    Dataset generate_awgn_dataset(const AwgnDatasetConfig& config);
    void save_dataset_csv(const Dataset& dataset,const std::filesystem::path& path);
    Dataset load_dataset_csv(const std::filesystem::path& path);
    void save_dataset_binary(const Dataset& dataset,const std::filesystem::path& path);
    Dataset load_dataset_binary(const std::filesystem::path& path);

    struct UnifiedDecodeResult
    {
        std::vector<DecoderPath> paths;
        DecoderStats stats;
        std::map<Rate1Engine,std::size_t> engine_counts;
        std::size_t tie_fallbacks{0};
    };

    UnifiedDecodeResult decode_variant(
    DecoderVariant variant,
    const std::vector<double>& llrs,
    const std::vector<bool>& frozen,
    std::size_t list_size,
    const PolarAlgorithmIR& portfolio_algorithm={});

    bool paths_equivalent(
    const std::vector<DecoderPath>& lhs,
    const std::vector<DecoderPath>& rhs,
    double tolerance=1e-10);

    struct ExperimentMetrics
    {
        DecoderVariant variant{DecoderVariant::SCL};
        std::size_t length{0};
        std::size_t information_bits{0};
        std::size_t payload_bits{0};
        std::size_t crc_bits{0};
        std::size_t list_size{0};
        double ebn0_db{0.0};
        std::size_t frames{0};
        std::size_t frame_errors{0};
        std::size_t bit_errors{0};
        std::size_t compared_frames{0};
        std::size_t list_mismatches{0};
        double mean_us{0.0};
        double p50_us{0.0};
        double p95_us{0.0};
        double p99_us{0.0};
        double max_us{0.0};
        double throughput_frames_per_second{0.0};
        double mean_leaf_steps{0.0};
        double mean_path_splits{0.0};
        double mean_prune_events{0.0};
        double mean_rate0_nodes{0.0};
        double mean_rate1_nodes{0.0};
        double mean_repetition_nodes{0.0};
        double mean_spc_nodes{0.0};
        double mean_special_candidates{0.0};
        double mean_rate1_candidates{0.0};
        double mean_frontier_expansions{0.0};
        double mean_rate1_fallbacks{0.0};
        double mean_tie_fallbacks{0.0};
        double mean_exhaustive_nodes{0.0};
        double mean_lazy_exact_nodes{0.0};
        double mean_fast_sscl_nodes{0.0};
        double mean_certified_frontier_nodes{0.0};

        double fer() const;
        double ber() const;
        double mismatch_rate() const;
    };

    struct ExperimentRunConfig
    {
        std::size_t list_size{8};
        bool compare_with_scl{true};
        std::size_t warmup_frames{8};
        std::size_t timing_repeats{1};
    };

    ExperimentMetrics run_decoder_experiment(
    const Dataset& dataset,
    DecoderVariant variant,
    const ExperimentRunConfig& run_config={},
    const PolarAlgorithmIR& portfolio_algorithm={});

    void write_metrics_csv(
    const std::vector<ExperimentMetrics>& metrics,
    const std::filesystem::path& path);
}
// namespace polar_scl
