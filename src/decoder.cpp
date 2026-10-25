#include "polar_scl/decoder.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace polar_scl
{
    namespace
    {
        double f_update(double left, double right)
        {
            const double sign = left * right >= 0.0 ? 1.0 : -1.0;
            return sign * std::min(std::abs(left), std::abs(right));
        }

        double g_update(double left, double right, Bit partial_sum)
        {
            return right + (1.0 - 2.0 * static_cast<double>(partial_sum)) * left;
        }

        double metric_increment(double llr, Bit bit)
        {
            const Bit hard = llr < 0.0 ? Bit{1} : Bit{0};
            return bit == hard ? 0.0 : std::abs(llr);
        }

        bool is_power_of_two(std::size_t x)
        {
            return x != 0 && (x & (x - 1)) == 0;
        }
    }

    // namespace

    bool decoder_path_less(const DecoderPath& lhs, const DecoderPath& rhs)
    {
        if(lhs.metric < rhs.metric) return true;
        if(rhs.metric < lhs.metric) return false;
        return std::lexicographical_compare(lhs.u.begin(), lhs.u.end(), rhs.u.begin(), rhs.u.end());
    }

    void DecoderStats::include_rate1(const SearchStats& s)
    {
        ++rate1_nodes;
        rate1_candidates += s.candidates_materialized;
        rate1_frontier_expansions += s.frontier_expansions;
        rate1_fallbacks += s.fallback_used ? 1 : 0;
    }

    std::vector<Bit> polar_transform(const std::vector<Bit>& bits)
    {
        if(!is_power_of_two(bits.size())) throw std::invalid_argument("polar_transform length must be a positive power of two");
        std::vector<Bit> values = bits;
        for(Bit bit : values) if(bit > 1) throw std::invalid_argument("polar_transform accepts binary values only");
        for(std::size_t step = 1; step < values.size(); step *= 2)
        {
            const std::size_t block = step * 2;
            for(std::size_t start = 0; start < values.size(); start += block)
            for(std::size_t offset = 0; offset < step; ++offset)
            values[start + offset] ^= values[start + step + offset];
        }

        return values;
    }

    std::vector<double> node_alpha(
    const std::vector<double>& channel_llrs,
    const std::vector<Bit>& prefix,
    std::size_t target_start,
    std::size_t target_size)
    {
        const std::size_t n = channel_llrs.size();
        if(!is_power_of_two(n)) throw std::invalid_argument("channel LLR length must be a positive power of two");
        if(!is_power_of_two(target_size)) throw std::invalid_argument("target_size must be a positive power of two");
        if(target_start + target_size > n) throw std::invalid_argument("target node lies outside the decoding tree");
        if(target_start % target_size != 0) throw std::invalid_argument("target node is not aligned to its size");
        if(prefix.size() != target_start) throw std::invalid_argument("prefix must contain decisions before target_start");

        std::function<std::vector<double>(const std::vector<double>&, std::size_t)> descend;
        descend = [&](const std::vector<double>& alpha, std::size_t start) -> std::vector<double>
        {
            const std::size_t size = alpha.size();
            if(start == target_start && size == target_size) return alpha;
            const std::size_t half = size / 2;
            const std::size_t middle = start + half;
            std::vector<double> left(alpha.begin(), alpha.begin() + static_cast<std::ptrdiff_t>(half));
            std::vector<double> right(alpha.begin() + static_cast<std::ptrdiff_t>(half), alpha.end());
            std::vector<double> child(half);
            if(target_start < middle)
            {
                for(std::size_t i = 0; i < half; ++i) child[i] = f_update(left[i], right[i]);
                return descend(child, start);
            }

            if(middle > prefix.size()) throw std::invalid_argument("right-child LLR requested before left sibling was decoded");
            std::vector<Bit> left_u(prefix.begin() + static_cast<std::ptrdiff_t>(start), prefix.begin() + static_cast<std::ptrdiff_t>(middle));
            auto beta_left = polar_transform(left_u);
            for(std::size_t i = 0; i < half; ++i) child[i] = g_update(left[i], right[i], beta_left[i]);
            return descend(child, middle);
        };
        return descend(channel_llrs, 0);
    }

    std::vector<DecoderPath> prune_paths(std::vector<DecoderPath> paths,std::size_t list_size)
    {
        if(paths.size()>list_size)
        {
            auto middle=paths.begin()+static_cast<std::ptrdiff_t>(list_size);
            std::nth_element(paths.begin(),middle,paths.end(),decoder_path_less);
            paths.resize(list_size);
        }
        std::sort(paths.begin(),paths.end(),decoder_path_less);
        return paths;
    }

    std::vector<DecoderPath> decode_leaf_range(
    const std::vector<double>& llrs,
    std::vector<DecoderPath> paths,
    const std::vector<bool>& frozen_flags,
    std::size_t start,
    std::size_t size,
    std::size_t list_size,
    DecoderStats& stats)
    {
        for(std::size_t index = start; index < start + size; ++index)
        {
            const bool is_frozen = frozen_flags[index];
            std::vector<DecoderPath> expanded;
            expanded.reserve(paths.size() * (is_frozen ? 1 : 2));
            for(const auto& path : paths)
            {
                const double llr = node_alpha(llrs, path.u, index, 1).front();
                const int count = is_frozen ? 1 : 2;
                for(int b = 0; b < count; ++b)
                {
                    const Bit bit = static_cast<Bit>(b);
                    auto u = path.u;
                    u.push_back(bit);
                    expanded.push_back({path.metric + metric_increment(llr, bit), std::move(u)});
                }
            }

            ++stats.leaf_steps;
            if(!is_frozen)
            {
                stats.path_splits += paths.size();
                ++stats.sorts;
            }

            paths = prune_paths(std::move(expanded), list_size);
        }

        return paths;
    }

    void validate_decoder_inputs(const std::vector<double>& llrs, const std::vector<bool>& frozen, std::size_t list_size)
    {
        if(!is_power_of_two(llrs.size())) throw std::invalid_argument("decoder length must be a positive power of two");
        if(llrs.size() != frozen.size()) throw std::invalid_argument("channel LLRs and frozen pattern must have equal length");
        if(list_size < 1) throw std::invalid_argument("list_size must be positive");
        for(double value : llrs) if(!std::isfinite(value)) throw std::invalid_argument("all channel LLRs must be finite");
    }

    DecodeResult scl_decode(const std::vector<double>& channel_llrs, const std::vector<bool>& frozen, std::size_t list_size)
    {
        validate_decoder_inputs(channel_llrs, frozen, list_size);
        DecoderStats stats;
        std::vector<DecoderPath> paths{{0.0, {}}};
        paths = decode_leaf_range(channel_llrs, std::move(paths), frozen, 0, channel_llrs.size(), list_size, stats);
        return {std::move(paths), stats};
    }

    std::vector<DecodeOp> rate1_schedule(const std::vector<bool>& frozen)
    {
        if(!is_power_of_two(frozen.size())) throw std::invalid_argument("frozen pattern length must be a positive power of two");
        std::vector<DecodeOp> operations;
        std::function<void(std::size_t, std::size_t)> visit;
        visit = [&](std::size_t start, std::size_t size)
        {
            bool all_frozen = true, any_frozen = false;
            for(std::size_t i = start; i < start + size; ++i)
            {
                all_frozen &= frozen[i];
                any_frozen |= frozen[i];
            }

            if(size > 1 && all_frozen)
            {
                operations.push_back({DecodeOp::Kind::Rate0, start, size});
                return;
            }

            if(size > 1 && !any_frozen && size <= 64)
            {
                operations.push_back({DecodeOp::Kind::Rate1, start, size});
                return;
            }

            if(size == 1)
            {
                operations.push_back({frozen[start] ? DecodeOp::Kind::Frozen : DecodeOp::Kind::Information, start, 1});
                return;
            }

            const auto half = size / 2;
            visit(start, half);
            visit(start + half, half);
        };
        visit(0, frozen.size());
        return operations;
    }

    std::vector<DecodeOp> enhanced_schedule(const std::vector<bool>& frozen)
    {
        if(!is_power_of_two(frozen.size()))
        {
            throw std::invalid_argument("frozen pattern length must be a positive power of two");
        }
        std::vector<DecodeOp> operations;
        std::function<void(std::size_t,std::size_t)> visit;
        visit=[&](std::size_t start,std::size_t size)
        {
            bool all_frozen=true;
            bool any_frozen=false;
            for(std::size_t index=start;index<start+size;++index)
            {
                all_frozen&=frozen[index];
                any_frozen|=frozen[index];
            }
            if(size>1&&all_frozen)
            {
                operations.push_back({DecodeOp::Kind::Rate0,start,size});
                return;
            }
            if(size>1&&!any_frozen&&size<=64)
            {
                operations.push_back({DecodeOp::Kind::Rate1,start,size});
                return;
            }
            if(size>1&&size<=64)
            {
                bool repetition=true;
                for(std::size_t index=start;index<start+size-1;++index)
                {
                    repetition&=frozen[index];
                }
                repetition&=!frozen[start+size-1];
                if(repetition)
                {
                    operations.push_back({DecodeOp::Kind::Repetition,start,size});
                    return;
                }
                bool spc=frozen[start];
                for(std::size_t index=start+1;index<start+size;++index)
                {
                    spc&=!frozen[index];
                }
                if(spc)
                {
                    operations.push_back({DecodeOp::Kind::SPC,start,size});
                    return;
                }
            }
            if(size==1)
            {
                operations.push_back({frozen[start]?DecodeOp::Kind::Frozen:DecodeOp::Kind::Information,start,1});
                return;
            }
            const std::size_t half=size/2;
            visit(start,half);
            visit(start+half,half);
        };
        visit(0,frozen.size());
        return operations;
    }

    DecodeResult rate1_scl_decode(
    const std::vector<double>& channel_llrs,
    const std::vector<bool>& frozen,
    std::size_t list_size,
    const Rate1Policy& policy)
    {
        validate_decoder_inputs(channel_llrs, frozen, list_size);
        std::vector<DecoderPath> paths{{0.0, {}}};
        DecoderStats stats;

        for(const auto& op : rate1_schedule(frozen))
        {
            for(const auto& path : paths) if(path.u.size() != op.start) throw std::runtime_error("decode schedule is not left-to-right");

            if(op.kind == DecodeOp::Kind::Rate0)
            {
                ++stats.rate0_nodes;
                std::vector<DecoderPath> updated;
                for(const auto& path : paths)
                {
                    const auto alpha = node_alpha(channel_llrs, path.u, op.start, op.size);
                    double penalty = 0.0;
                    for(double value : alpha) if(value < 0.0) penalty += std::abs(value);
                    auto u = path.u;
                    u.insert(u.end(), op.size, Bit{0}

                    );
                    updated.push_back({path.metric + penalty, std::move(u)});
                }

                paths = prune_paths(std::move(updated), list_size);
                continue;
            }

            if(op.kind == DecodeOp::Kind::Rate1)
            {
                std::vector<ParentPath> parents;
                parents.reserve(paths.size());
                for(std::size_t i = 0; i < paths.size(); ++i)
                parents.push_back({i, paths[i].metric, node_alpha(channel_llrs, paths[i].u, op.start, op.size)});

                if(rate1_metric_tie_risk(parents))
                {
                    paths=decode_leaf_range(channel_llrs,std::move(paths),frozen,op.start,op.size,list_size,stats);
                    continue;
                }

                CandidateRankFn rank_fn = [&](const Rate1Candidate& candidate)
                {
                    auto local_u = polar_transform(candidate.bits);
                    std::vector<std::uint64_t> tie;
                    tie.reserve(paths[candidate.parent_id].u.size() + local_u.size());
                    for(Bit b : paths[candidate.parent_id].u) tie.push_back(b);
                    for(Bit b : local_u) tie.push_back(b);
                    return GeneralRankKey{candidate.path_metric, std::move(tie)};
                };
                HiddenRankFn hidden_fn = [&](const ParentPath& parent, double metric)
                {
                    std::vector<std::uint64_t> tie;
                    tie.reserve(paths[parent.parent_id].u.size() + op.size);
                    for(Bit b : paths[parent.parent_id].u) tie.push_back(b);
                    tie.insert(tie.end(), op.size, 0);
                    return GeneralRankKey{metric, std::move(tie)};
                };
                auto node = adaptive_rate1_top_l(parents, list_size, policy, rank_fn, hidden_fn);
                stats.include_rate1(node.stats);
                if(node.stats.fallback_used)
                {
                    paths = decode_leaf_range(channel_llrs, std::move(paths), frozen, op.start, op.size, list_size, stats);
                    continue;
                }

                std::vector<DecoderPath> expanded;
                for(const auto& candidate : node.candidates)
                {
                    auto u = paths[candidate.parent_id].u;
                    auto local_u = polar_transform(candidate.bits);
                    u.insert(u.end(), local_u.begin(), local_u.end());
                    expanded.push_back({candidate.path_metric, std::move(u)});
                }

                paths = prune_paths(std::move(expanded), list_size);
                ++stats.sorts;
                continue;
            }

            paths = decode_leaf_range(channel_llrs, std::move(paths), frozen, op.start, 1, list_size, stats);
        }

        return {std::move(paths), stats};
    }
}

// namespace polar_scl
