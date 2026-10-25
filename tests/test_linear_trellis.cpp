#include "test_util.hpp"

#include <random>

using namespace polar_scl;

namespace
{
    std::vector<Rate1Candidate> exhaustive_linear(
        const std::vector<ParentPath>& parents, std::size_t list_size,
        const LinearConstraintIR& constraint)
    {
        std::vector<Rate1Candidate> all;
        for(const auto& parent : parents)
        {
            const auto limit = std::uint64_t{1} << parent.llrs.size();
            for(std::uint64_t mask = 0; mask < limit; ++mask)
            {
                auto candidate = Rate1Candidate{};
                candidate.parent_id = parent.parent_id;
                candidate.flip_mask = mask;
                candidate.bits.resize(parent.llrs.size());
                for(std::size_t i = 0; i < parent.llrs.size(); ++i)
                {
                    const auto hard = static_cast<Bit>(parent.llrs[i] < 0.0);
                    const auto flip = static_cast<Bit>((mask >> i) & 1ULL);
                    candidate.bits[i] = static_cast<Bit>(hard ^ flip);
                    if(flip) candidate.penalty += std::abs(parent.llrs[i]);
                }
                bool valid = true;
                for(auto row : constraint.parity_checks)
                {
                    unsigned parity = 0;
                    for(std::size_t i = 0; i < parent.llrs.size(); ++i)
                        if(((row >> i) & 1ULL) != 0ULL) parity ^= candidate.bits[i];
                    valid = valid && parity == 0;
                }
                if(valid)
                {
                    candidate.path_metric = parent.path_metric + candidate.penalty;
                    all.push_back(std::move(candidate));
                }
            }
        }
        std::sort(all.begin(), all.end(), [](const auto& a, const auto& b)
        { return local_rank_key(a) < local_rank_key(b); });
        if(all.size() > list_size) all.resize(list_size);
        return all;
    }
}

int main()
{
    return run_test_main("linear_trellis", []
    {
        const LinearConstraintIR redundant{"overlap", 8, {0x55, 0xAA, 0xFF}};
        CHECK_EQ(redundant.rank(), 2u);
        CHECK_EQ(redundant.state_count(), 4u);

        TrellisSearchReport report;
        const auto searched = search_trellis_ir(
            LinearConstraintIR{"interleaved", 8, {0x55, 0xAA}}, 256, 73, &report);
        CHECK_EQ(report.unique_structures, 255u);
        CHECK_TRUE(report.selected_score < report.baseline_score);
        CHECK_EQ(report.baseline_score, (std::tuple<std::size_t,std::size_t,std::size_t>{26,4,31}));
        CHECK_EQ(report.selected_score, (std::tuple<std::size_t,std::size_t,std::size_t>{15,2,25}));
        auto selected_order = report.selected_order;
        std::sort(selected_order.begin(), selected_order.end());
        CHECK_EQ(selected_order, (std::vector<std::size_t>{0,1,2,3,4,5,6,7}));

        std::mt19937_64 rng(701);
        std::uniform_int_distribution<int> row_dist(1, 255), llr_dist(-5, 5);
        for(int matrix_index = 0; matrix_index < 12; ++matrix_index)
        {
            LinearConstraintIR constraint{"random", 8,
                {static_cast<std::uint64_t>(row_dist(rng)),
                 static_cast<std::uint64_t>(row_dist(rng)),
                 static_cast<std::uint64_t>(row_dist(rng))}};
            auto trellis = search_trellis_ir(constraint, 32, 1000u + matrix_index);
            for(std::size_t list_size : {1u, 2u, 4u, 8u})
            {
                for(int sample = 0; sample < 5; ++sample)
                {
                    std::vector<ParentPath> parents;
                    for(std::size_t p = 0; p < 2; ++p)
                    {
                        std::vector<double> llrs(8);
                        for(auto& value : llrs) value = llr_dist(rng);
                        parents.push_back({p, static_cast<double>(p), std::move(llrs)});
                    }
                    const auto expected = exhaustive_linear(parents, list_size, constraint);
                    const auto actual = linear_trellis_top_l(parents, list_size, trellis);
                    check_candidates_equal(actual.candidates, expected);
                }
            }
        }

        CHECK_THROW(search_trellis_ir(redundant, 1, 0, nullptr, 2));
        CHECK_THROW((TrellisIR{redundant, {0,0,0,0,0,0,0,0}}).validate());
    });
}
