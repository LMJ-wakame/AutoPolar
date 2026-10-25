#pragma once

#include "polar_scl/ir.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace polar_scl
{
    using Bit = std::uint8_t;

    struct ParentPath
    {
        std::size_t parent_id{};
        double path_metric{};
        std::vector<double> llrs;
    };

    struct Rate1Candidate
    {
        std::size_t parent_id{};
        double path_metric{};
        double penalty{};
        std::uint64_t flip_mask{};
        std::vector<Bit> bits;
    };

    struct LocalRankKey
    {
        double metric{};
        std::size_t parent_id{};
        std::uint64_t flip_mask{};

        auto operator<=>(const LocalRankKey&) const = default;
    };

    struct GeneralRankKey
    {
        double metric{};
        std::vector<std::uint64_t> tie;
    };

    bool operator<(const GeneralRankKey& lhs, const GeneralRankKey& rhs);
    bool operator>(const GeneralRankKey& lhs, const GeneralRankKey& rhs);
    bool operator<=(const GeneralRankKey& lhs, const GeneralRankKey& rhs);
    bool operator==(const GeneralRankKey& lhs, const GeneralRankKey& rhs);

    struct SearchStats
    {
        std::string algorithm;
        std::size_t candidates_materialized{0};
        std::size_t subset_states_pushed{0};
        std::size_t subset_states_popped{0};
        std::size_t global_heap_pushes{0};
        std::size_t global_heap_pops{0};
        std::size_t frontier_expansions{0};
        std::size_t certificate_checks{0};
        bool fallback_used{false};
        std::string fallback_reason;
        std::optional<GeneralRankKey> certificate_worst_key;
        std::optional<GeneralRankKey> certificate_hidden_key;
        std::map<std::size_t, std::size_t> per_parent_active_bits;
        std::size_t parity_rejections{}, syndrome_rejections{};
        std::size_t metric_evaluations{}, constraint_checks{}, feasible_candidates{}, dp_state_relaxations{};
        std::size_t primitive_work() const
        {
            return metric_evaluations + constraint_checks + subset_states_pushed + subset_states_popped
                + global_heap_pushes + global_heap_pops + dp_state_relaxations;
        }
    };

    struct Rate1SearchResult
    {
        std::vector<Rate1Candidate> candidates;
        SearchStats stats;
    };

    using CandidateRankFn = std::function<GeneralRankKey(const Rate1Candidate&)>;
    using HiddenRankFn = std::function<GeneralRankKey(const ParentPath&, double)>;

    LocalRankKey local_rank_key(const Rate1Candidate& candidate);
    GeneralRankKey general_local_rank_key(const Rate1Candidate& candidate);

    void validate_parent_paths(const std::vector<ParentPath>& parents, std::size_t list_size);
    bool rate1_metric_tie_risk(const std::vector<ParentPath>& parents,double tolerance=1e-9);

    Rate1SearchResult exhaustive_rate1_top_l(
    const std::vector<ParentPath>& parents,
    std::size_t list_size);

    Rate1SearchResult lazy_exact_rate1_top_l(
    const std::vector<ParentPath>& parents,
    std::size_t list_size);

    Rate1SearchResult fast_sscl_bounded_rate1_top_l(
    const std::vector<ParentPath>& parents,
    std::size_t list_size);

    Rate1SearchResult adaptive_rate1_top_l(
    const std::vector<ParentPath>& parents,
    std::size_t list_size,
    const Rate1Policy& policy = {},
    CandidateRankFn rank_key = {},
    HiddenRankFn hidden_rank_key = {});
}

// namespace polar_scl
