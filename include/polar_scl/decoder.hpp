#pragma once

#include "polar_scl/ir.hpp"
#include "polar_scl/rate1.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace polar_scl
{
    struct DecoderPath
    {
        double metric{};
        std::vector<Bit> u;
    };

    bool decoder_path_less(const DecoderPath& lhs, const DecoderPath& rhs);

    struct DecodeOp
    {
        enum class Kind
        {
            Frozen, Information, Rate0, Rate1, Repetition, SPC
        };
        Kind kind{Kind::Frozen};
        std::size_t start{};
        std::size_t size{};
    };

    struct DecoderStats
    {
        std::size_t leaf_steps{0};
        std::size_t path_splits{0};
        std::size_t sorts{0};
        std::size_t rate0_nodes{0};
        std::size_t rate1_nodes{0};
        std::size_t repetition_nodes{0};
        std::size_t spc_nodes{0};
        std::size_t special_candidates{0};
        std::size_t rate1_candidates{0};
        std::size_t rate1_frontier_expansions{0};
        std::size_t rate1_fallbacks{0};

        void include_rate1(const SearchStats& stats);
    };

    struct DecodeResult
    {
        std::vector<DecoderPath> paths;
        DecoderStats stats;
    };

    std::vector<Bit> polar_transform(const std::vector<Bit>& bits);

    std::vector<DecodeOp> rate1_schedule(const std::vector<bool>& frozen);
    std::vector<DecodeOp> enhanced_schedule(const std::vector<bool>& frozen);

    DecodeResult scl_decode(
    const std::vector<double>& channel_llrs,
    const std::vector<bool>& frozen,
    std::size_t list_size);

    DecodeResult rate1_scl_decode(
    const std::vector<double>& channel_llrs,
    const std::vector<bool>& frozen,
    std::size_t list_size,
    const Rate1Policy& policy = {});

    // Public helper used by the portfolio decoder.
    std::vector<double> node_alpha(
    const std::vector<double>& channel_llrs,
    const std::vector<Bit>& prefix,
    std::size_t target_start,
    std::size_t target_size);

    std::vector<DecoderPath> prune_paths(
    std::vector<DecoderPath> paths,
    std::size_t list_size);

    std::vector<DecoderPath> decode_leaf_range(
    const std::vector<double>& llrs,
    std::vector<DecoderPath> paths,
    const std::vector<bool>& frozen_flags,
    std::size_t start,
    std::size_t size,
    std::size_t list_size,
    DecoderStats& stats);

    void validate_decoder_inputs(
    const std::vector<double>& llrs,
    const std::vector<bool>& frozen,
    std::size_t list_size);
}

// namespace polar_scl
