#include "polar_scl/linear_trellis.hpp"
#include "polar_scl_lab/linear.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>

namespace polar_scl
{
    namespace
    {
        unsigned highest_bit(std::uint64_t value)
        {
            return 63u - static_cast<unsigned>(std::countl_zero(value));
        }

        std::vector<std::uint64_t> compressed_columns(const LinearConstraintIR& constraint)
        {
            constraint.validate();
            std::vector<std::uint64_t> raw(constraint.num_variables, 0);
            for(std::size_t row = 0; row < constraint.parity_checks.size(); ++row)
            {
                for(std::size_t col = 0; col < constraint.num_variables; ++col)
                {
                    if(((constraint.parity_checks[row] >> col) & 1ULL) != 0ULL)
                        raw[col] |= (1ULL << row);
                }
            }

            std::vector<std::uint64_t> echelon(64, 0);
            for(auto value : raw)
            {
                while(value != 0)
                {
                    const auto pivot = highest_bit(value);
                    if(echelon[pivot] == 0)
                    {
                        echelon[pivot] = value;
                        break;
                    }
                    value ^= echelon[pivot];
                }
            }

            std::vector<std::pair<unsigned, std::uint64_t>> basis;
            for(unsigned pivot = 0; pivot < 64; ++pivot)
                if(echelon[pivot] != 0) basis.emplace_back(pivot, echelon[pivot]);
            std::reverse(basis.begin(), basis.end());

            std::vector<std::uint64_t> result;
            result.reserve(raw.size());
            for(auto value : raw)
            {
                std::uint64_t coordinates = 0;
                for(std::size_t index = 0; index < basis.size(); ++index)
                {
                    const auto [pivot, vector] = basis[index];
                    if(((value >> pivot) & 1ULL) != 0ULL)
                    {
                        value ^= vector;
                        coordinates |= (1ULL << index);
                    }
                }
                if(value != 0) throw std::logic_error("failed to compile constraint coordinates");
                result.push_back(coordinates);
            }
            return result;
        }

        using StateSet = std::set<std::uint64_t>;

        std::size_t raw_column_rank(const LinearConstraintIR& constraint)
        {
            std::vector<std::uint64_t> echelon(64, 0);
            std::size_t count = 0;
            for(std::size_t col = 0; col < constraint.num_variables; ++col)
            {
                std::uint64_t value = 0;
                for(std::size_t row = 0; row < constraint.parity_checks.size(); ++row)
                    if(((constraint.parity_checks[row] >> col) & 1ULL) != 0ULL)
                        value |= (1ULL << row);
                while(value != 0)
                {
                    const auto pivot = highest_bit(value);
                    if(echelon[pivot] == 0)
                    {
                        echelon[pivot] = value;
                        ++count;
                        break;
                    }
                    value ^= echelon[pivot];
                }
            }
            return count;
        }

        std::vector<StateSet> suffix_spaces(const TrellisIR& trellis, std::size_t max_states)
        {
            const auto cols = trellis.constraint.columns();
            const auto n = trellis.order.size();
            std::vector<StateSet> suffix(n + 1);
            suffix[n].insert(0);
            for(std::size_t depth = n; depth-- > 0;)
            {
                suffix[depth] = suffix[depth + 1];
                const auto column = cols[trellis.order[depth]];
                std::vector<std::uint64_t> additions;
                additions.reserve(suffix[depth].size());
                for(auto state : suffix[depth]) additions.push_back(state ^ column);
                suffix[depth].insert(additions.begin(), additions.end());
                if(suffix[depth].size() > max_states)
                    throw std::length_error("trellis suffix-state limit exceeded");
            }
            return suffix;
        }

        std::vector<std::uint64_t> hard_syndrome(const ParentPath& parent,
                                                  const LinearConstraintIR& constraint)
        {
            const auto cols = constraint.columns();
            std::uint64_t result = 0;
            for(std::size_t i = 0; i < parent.llrs.size(); ++i)
                if(parent.llrs[i] < 0.0) result ^= cols[i];
            return {result};
        }

        Rate1Candidate candidate_for(const ParentPath& parent, std::uint64_t mask)
        {
            return polar_scl_lab::make_candidate(parent, mask);
        }

        bool candidate_less(const Rate1Candidate& a, const Rate1Candidate& b)
        {
            return local_rank_key(a) < local_rank_key(b);
        }
    }

    void LinearConstraintIR::validate() const
    {
        if(name.find_first_not_of(" \t\r\n")==std::string::npos) throw std::invalid_argument("constraint name must not be empty");
        if(num_variables == 0 || num_variables > 64)
            throw std::invalid_argument("constraint supports 1..64 variables");
        if(parity_checks.empty() || parity_checks.size() > 64)
            throw std::invalid_argument("constraint supports 1..64 parity checks");
        const auto valid_mask = num_variables == 64 ? ~0ULL : ((1ULL << num_variables) - 1ULL);
        for(auto row : parity_checks)
            if(row == 0 || (row & ~valid_mask) != 0)
                throw std::invalid_argument("parity-check row is zero or out of range");
    }

    std::size_t LinearConstraintIR::rank() const
    {
        validate();
        return raw_column_rank(*this);
    }

    std::vector<std::uint64_t> LinearConstraintIR::raw_columns() const
    {
        validate();std::vector<std::uint64_t> raw(num_variables);
        for(std::size_t row=0;row<parity_checks.size();++row)
            for(std::size_t col=0;col<num_variables;++col)
                if((parity_checks[row]>>col)&1ULL)raw[col]|=1ULL<<row;
        return raw;
    }

    std::size_t LinearConstraintIR::state_count() const
    {
        const auto r = rank();
        if(r >= std::numeric_limits<std::size_t>::digits)
            throw std::length_error("compiled syndrome space is not addressable");
        return std::size_t{1} << r;
    }

    std::vector<std::uint64_t> LinearConstraintIR::columns() const
    {
        return compressed_columns(*this);
    }

    std::uint64_t LinearConstraintIR::syndrome(std::uint64_t word_mask) const
    {
        if(num_variables < 64 && (word_mask >> num_variables) != 0)
            throw std::invalid_argument("word mask is outside constraint domain");
        const auto cols = columns();
        std::uint64_t result = 0;
        for(std::size_t i = 0; i < num_variables; ++i)
            if(((word_mask >> i) & 1ULL) != 0ULL) result ^= cols[i];
        return result;
    }

    void TrellisIR::validate() const
    {
        constraint.validate();
        if(order.size() != constraint.num_variables)
            throw std::invalid_argument("trellis order has incorrect length");
        auto sorted = order;
        std::sort(sorted.begin(), sorted.end());
        for(std::size_t i = 0; i < sorted.size(); ++i)
            if(sorted[i] != i) throw std::invalid_argument("trellis order must be a permutation");
    }

    std::vector<std::size_t> TrellisIR::width_profile() const
    {
        validate();
        if(constraint.rank() > 20)
            throw std::length_error("trellis width analysis is limited to at most 2^20 states");
        const auto cols = constraint.columns();
        std::vector<StateSet> prefix(order.size() + 1), suffix(order.size() + 1);
        prefix[0].insert(0);
        suffix[order.size()].insert(0);
        for(std::size_t d = 0; d < order.size(); ++d)
        {
            prefix[d + 1] = prefix[d];
            const auto column = cols[order[d]];
            for(auto state : prefix[d]) prefix[d + 1].insert(state ^ column);
        }
        for(std::size_t d = order.size(); d-- > 0;)
        {
            suffix[d] = suffix[d + 1];
            const auto column = cols[order[d]];
            for(auto state : suffix[d + 1]) suffix[d].insert(state ^ column);
        }
        std::vector<std::size_t> widths(order.size() + 1);
        for(std::size_t d = 0; d <= order.size(); ++d)
        {
            for(auto state : prefix[d]) if(suffix[d].contains(state)) ++widths[d];
        }
        return widths;
    }

    auto TrellisIR::score() const -> std::tuple<std::size_t, std::size_t, std::size_t>
    {
        const auto widths = width_profile();
        const auto sum_width = std::accumulate(widths.begin(), widths.end(), std::size_t{0});
        const auto max_width = *std::max_element(widths.begin(), widths.end());
        const auto cols = constraint.columns();
        std::vector<StateSet> suffix(order.size() + 1);
        suffix.back().insert(0);
        for(std::size_t d = order.size(); d-- > 0;)
        {
            suffix[d] = suffix[d + 1];
            const auto column = cols[order[d]];
            for(auto state : suffix[d + 1]) suffix[d].insert(state ^ column);
        }
        std::size_t suffix_entries = 0;
        for(const auto& states : suffix) suffix_entries += states.size();
        return {sum_width, max_width, suffix_entries};
    }

    TrellisIR search_trellis_ir(const LinearConstraintIR& constraint, std::size_t trials,
                                std::uint64_t seed, TrellisSearchReport* report,
                                std::size_t max_states)
    {
        constraint.validate();
        if(trials == 0) throw std::invalid_argument("trials must be positive");
        if(max_states == 0 || constraint.rank() >= std::numeric_limits<std::size_t>::digits ||
           (std::size_t{1} << constraint.rank()) > max_states)
            throw std::length_error("constraint state space exceeds configured search limit");
        TrellisIR baseline;
        baseline.constraint = constraint;
        baseline.order.resize(constraint.num_variables);
        std::iota(baseline.order.begin(), baseline.order.end(), 0);
        std::set<std::vector<std::size_t>> orders{baseline.order};
        std::mt19937_64 rng(seed);
        for(std::size_t i = 1; i < trials; ++i)
        {
            auto order = baseline.order;
            std::shuffle(order.begin(), order.end(), rng);
            orders.insert(std::move(order));
        }
        auto best = baseline;
        auto best_score = best.score();
        for(const auto& order : orders)
        {
            TrellisIR candidate{constraint, order};
            const auto score = candidate.score();
            if(score < best_score || (score == best_score && candidate.order < best.order))
            {
                best = std::move(candidate);
                best_score = score;
            }
        }
        if(report)
        {
            report->seed = seed;
            report->trials = trials;
            report->unique_structures = orders.size();
            report->selected_order = best.order;
            report->baseline_score = baseline.score();
            report->selected_score = best_score;
        }
        return best;
    }

    Rate1SearchResult linear_trellis_top_l(const std::vector<ParentPath>& parents,
                                            std::size_t list_size, const TrellisIR& trellis,
                                            LinearTrellisStats* trellis_stats)
    {
        trellis.validate();
        validate_parent_paths(parents, list_size);
        for(const auto& parent : parents)
            if(parent.llrs.size() != trellis.constraint.num_variables)
                throw std::invalid_argument("parent length does not match constraint");

        constexpr std::size_t max_states = 1u << 20;
        const auto suffix = suffix_spaces(trellis, max_states);
        const auto columns = trellis.constraint.columns();
        Rate1SearchResult result;
        result.stats.algorithm = "linear_sparse_trellis_top_l";
        LinearTrellisStats local_stats;
        for(const auto& parent : parents)
        {
            const auto target = hard_syndrome(parent, trellis.constraint).front();
            std::map<std::uint64_t, std::vector<std::uint64_t>> states;
            states[0].push_back(0);
            for(std::size_t depth = 0; depth < trellis.order.size(); ++depth)
            {
                std::map<std::uint64_t, std::vector<std::uint64_t>> next;
                const auto variable = trellis.order[depth];
                for(const auto& [syndrome, masks] : states)
                {
                    for(auto mask : masks)
                    {
                        for(unsigned bit = 0; bit < 2; ++bit)
                        {
                            ++local_stats.transitions;
                            const auto dest = syndrome ^ (bit ? columns[variable] : 0ULL);
                            if(!suffix[depth + 1].contains(target ^ dest)) continue;
                            auto& bucket = next[dest];
                            bucket.push_back(bit ? (mask | (1ULL << variable)) : mask);
                        }
                    }
                }
                for(auto& [state, masks] : next)
                {
                    (void)state;
                    std::sort(masks.begin(), masks.end(), [&](auto a, auto b)
                    { return candidate_less(candidate_for(parent, a), candidate_for(parent, b)); });
                    if(masks.size() > list_size) masks.resize(list_size);
                }
                states = std::move(next);
            }
            const auto found = states.find(target);
            if(found != states.end())
                for(auto mask : found->second)
                    result.candidates.push_back(candidate_for(parent, mask));
        }
        std::sort(result.candidates.begin(), result.candidates.end(), candidate_less);
        if(result.candidates.size() > list_size) result.candidates.resize(list_size);
        result.stats.candidates_materialized = result.candidates.size();
        for(const auto& layer : suffix) local_stats.suffix_state_entries += layer.size();
        if(trellis_stats) *trellis_stats = local_stats;
        return result;
    }
}
