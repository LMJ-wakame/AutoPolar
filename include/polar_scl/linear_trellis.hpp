#pragma once

#include "polar_scl/rate1.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace polar_scl
{
    // Binary parity-check IR. Rows use bit i for variable i; at most 64
    // variables and checks are represented by the current portable backend.
    struct LinearConstraintIR
    {
        std::string name;
        std::size_t num_variables{};
        std::vector<std::uint64_t> parity_checks;

        std::size_t rank() const;
        std::vector<std::uint64_t> raw_columns() const;
        std::size_t state_count() const;
        std::vector<std::uint64_t> columns() const;
        std::uint64_t syndrome(std::uint64_t word_mask) const;
        void validate() const;
    };

    struct TrellisIR
    {
        LinearConstraintIR constraint;
        std::vector<std::size_t> order;

        void validate() const;
        std::vector<std::size_t> width_profile() const;
        auto score() const -> std::tuple<std::size_t, std::size_t, std::size_t>;
    };

    struct TrellisSearchReport
    {
        std::uint64_t seed{};
        std::size_t trials{};
        std::size_t unique_structures{};
        std::vector<std::size_t> selected_order;
        std::tuple<std::size_t, std::size_t, std::size_t> baseline_score{};
        std::tuple<std::size_t, std::size_t, std::size_t> selected_score{};
    };

    struct LinearTrellisStats
    {
        std::size_t transitions{};
        std::size_t suffix_state_entries{};
    };

    TrellisIR search_trellis_ir(
        const LinearConstraintIR& constraint,
        std::size_t trials,
        std::uint64_t seed,
        TrellisSearchReport* report = nullptr,
        std::size_t max_states = (1u << 20));

    Rate1SearchResult linear_trellis_top_l(
        const std::vector<ParentPath>& parents,
        std::size_t list_size,
        const TrellisIR& trellis,
        LinearTrellisStats* trellis_stats = nullptr);
}
