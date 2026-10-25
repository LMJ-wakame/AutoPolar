#pragma once

#include "polar_scl/algorithm_ir.hpp"
#include "polar_scl/decoder.hpp"
#include "polar_scl/rate1.hpp"
#include "polar_scl/selector.hpp"

#include <map>

namespace polar_scl
{
    struct PortfolioChoice
    {
        Rate1NodeIR ir;
        Rate1NodeState state;
    };

    PortfolioChoice choose_rate1_ir(
    const PolarAlgorithmIR& algorithm,
    const std::vector<ParentPath>& parents,
    std::size_t list_size);

    struct PortfolioRate1Result
    {
        std::vector<Rate1Candidate> candidates;
        SearchStats stats;
        Rate1NodeIR ir;
        Rate1NodeState state;
    };

    PortfolioRate1Result portfolio_rate1_top_l(
    const PolarAlgorithmIR& algorithm,
    const std::vector<ParentPath>& parents,
    std::size_t list_size);

    struct PortfolioDecoderStats
    {
        DecoderStats decoder;
        std::map<Rate1Engine, std::size_t> engine_counts;
        std::size_t tie_fallbacks{0};
    };

    struct PortfolioDecodeResult
    {
        std::vector<DecoderPath> paths;
        PortfolioDecoderStats stats;
    };

    struct PortfolioDecodeOptions
    {
        bool enable_rate1_tie_guard{true};
        bool enable_spc_cutoff_guard{true};
        bool enable_final_list_tie_guard{true};
    };

    PortfolioDecodeResult portfolio_scl_decode(
    const std::vector<double>& channel_llrs,
    const std::vector<bool>& frozen,
    std::size_t list_size,
    const PolarAlgorithmIR& algorithm={},
    const PortfolioDecodeOptions& options={});
}

// namespace polar_scl
