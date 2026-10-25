#include "polar_scl/portfolio.hpp"
#include "polar_scl/materialize.hpp"
#include "polar_scl/special_nodes.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace polar_scl
{
    namespace
    {
        bool metric_close(double lhs,double rhs,double tolerance=1e-9)
        {
            return std::abs(lhs-rhs)<=tolerance*std::max({1.0,std::abs(lhs),std::abs(rhs)});
        }

        bool has_close_path_metrics(const std::vector<DecoderPath>& paths)
        {
            for(std::size_t left=0;left<paths.size();++left)
            {
                for(std::size_t right=left+1;right<paths.size();++right)
                {
                    if(metric_close(paths[left].metric,paths[right].metric))
                    {
                        return true;
                    }
                }
            }
            return false;
        }
    }

    PortfolioChoice choose_rate1_ir(
        const PolarAlgorithmIR& algorithm,
        const std::vector<ParentPath>& parents,
        std::size_t list_size)
    {
        validate_algorithm_ir(algorithm);
        if(parents.empty())
        {
            throw std::invalid_argument("at least one Rate-1 parent is required");
        }

        Rate1NodeState basic;
        basic.node_size=parents.front().llrs.size();
        basic.list_size=list_size;
        basic.parent_count=parents.size();
        if(algorithm.selector_mode==SelectorMode::Fixed)
        {
            return {algorithm.fixed_rate1,basic};
        }

        const SelectorGenome& genome=algorithm.selector;
        if(genome.list_size_threshold>0&&list_size<=genome.list_size_threshold)
        {
            return {genome.ir_for_engine(genome.low_list_engine),basic};
        }

        if(basic.node_size<=genome.small_node_threshold)
        {
            return {genome.ir_for_engine(genome.small_engine),basic};
        }

        if(genome.easy_engine==genome.hard_engine)
        {
            return {genome.ir_for_engine(genome.easy_engine),basic};
        }

        Rate1NodeState state=build_rate1_state(parents,list_size);
        return {select_rate1_ir(genome,state),state};
    }

    PortfolioRate1Result portfolio_rate1_top_l(const PolarAlgorithmIR& algorithm, const std::vector<ParentPath>& parents, std::size_t list_size)
    {
        auto choice = choose_rate1_ir(algorithm, parents, list_size);
        auto result = execute_rate1_ir(choice.ir, parents, list_size);
        return {std::move(result.candidates), std::move(result.stats), choice.ir, choice.state};
    }

    PortfolioDecodeResult portfolio_scl_decode(
    const std::vector<double>& channel_llrs,
    const std::vector<bool>& frozen,
    std::size_t list_size,
    const PolarAlgorithmIR& algorithm,
    const PortfolioDecodeOptions& options)
    {
        validate_algorithm_ir(algorithm);
        validate_decoder_inputs(channel_llrs, frozen, list_size);
        std::vector<DecoderPath> paths{{0.0, {}}};
        PortfolioDecoderStats result_stats;
        auto& stats = result_stats.decoder;

        for(const auto& op:enhanced_schedule(frozen))
        {
            for(const auto& path : paths) if(path.u.size() != op.start) throw std::runtime_error("decode schedule is not left-to-right");
            if(op.kind == DecodeOp::Kind::Rate0)
            {
                ++stats.rate0_nodes;
                if(paths.size()>1)
                {
                    ++result_stats.tie_fallbacks;
                    paths=decode_leaf_range(
                        channel_llrs,
                        std::move(paths),
                        frozen,
                        op.start,
                        op.size,
                        list_size,
                        stats);
                    continue;
                }
                std::vector<DecoderPath> updated;
                for(const auto& path : paths)
                {
                    auto alpha = node_alpha(channel_llrs, path.u, op.start, op.size);
                    double penalty = 0.0;
                    for(double v : alpha) if(v < 0.0) penalty += std::abs(v);
                    auto u = path.u;
                    u.insert(u.end(), op.size, Bit{0}

                    );
                    updated.push_back({path.metric + penalty, std::move(u)});
                }

                if(options.enable_rate1_tie_guard&&has_close_path_metrics(updated))
                {
                    ++result_stats.tie_fallbacks;
                    paths=decode_leaf_range(
                        channel_llrs,
                        std::move(paths),
                        frozen,
                        op.start,
                        op.size,
                        list_size,
                        stats);
                    continue;
                }

                paths=prune_paths(std::move(updated),list_size);
                continue;
            }

            if(op.kind==DecodeOp::Kind::Repetition||op.kind==DecodeOp::Kind::SPC)
            {
                const bool enabled=algorithm.selector_mode==SelectorMode::Rule&&
                    (op.kind==DecodeOp::Kind::Repetition
                        ?algorithm.selector.enable_repetition
                        :(algorithm.selector.enable_spc&&list_size<=8));
                if(!enabled)
                {
                    paths=decode_leaf_range(channel_llrs,std::move(paths),frozen,op.start,op.size,list_size,stats);
                    continue;
                }
                std::vector<ParentPath> parents;
                parents.reserve(paths.size());
                for(std::size_t index=0;index<paths.size();++index)
                {
                    parents.push_back({index,paths[index].metric,node_alpha(channel_llrs,paths[index].u,op.start,op.size)});
                }
                if(options.enable_rate1_tie_guard&&rate1_metric_tie_risk(parents))
                {
                    ++result_stats.tie_fallbacks;
                    paths=decode_leaf_range(
                        channel_llrs,
                        std::move(paths),
                        frozen,
                        op.start,
                        op.size,
                        list_size,
                        stats);
                    continue;
                }
                CandidateRankFn rank_fn=[&](const Rate1Candidate& candidate)
                {
                    const auto local_u=polar_transform(candidate.bits);
                    std::vector<std::uint64_t> tie;
                    tie.reserve(paths[candidate.parent_id].u.size()+local_u.size());
                    for(Bit bit:paths[candidate.parent_id].u)
                    {
                        tie.push_back(bit);
                    }
                    for(Bit bit:local_u)
                    {
                        tie.push_back(bit);
                    }
                    return GeneralRankKey{candidate.path_metric,std::move(tie)};
                };
                Rate1SearchResult node;
                if(op.kind==DecodeOp::Kind::Repetition)
                {
                    node=repetition_node_top_l(parents,list_size+1,rank_fn);
                    ++stats.repetition_nodes;
                }
                else
                {
                    node=spc_node_top_l(parents,list_size+1,rank_fn);
                    ++stats.spc_nodes;
                }
                if(options.enable_spc_cutoff_guard)
                {
                    bool internal_tie=false;
                    for(std::size_t left=0;left<node.candidates.size()&&!internal_tie;++left)
                    {
                        for(std::size_t right=left+1;right<node.candidates.size();++right)
                        {
                            if(metric_close(
                                node.candidates[left].path_metric,
                                node.candidates[right].path_metric))
                            {
                                internal_tie=true;
                                break;
                            }
                        }
                    }
                    if(internal_tie)
                    {
                        ++result_stats.tie_fallbacks;
                        paths=decode_leaf_range(
                            channel_llrs,
                            std::move(paths),
                            frozen,
                            op.start,
                            op.size,
                            list_size,
                            stats);
                        continue;
                    }
                }
                if(node.candidates.size()>list_size)
                {
                    const double cutoff=node.candidates[list_size-1].path_metric;
                    const double next=node.candidates[list_size].path_metric;
                    if(options.enable_spc_cutoff_guard&&metric_close(cutoff,next))
                    {
                        ++result_stats.tie_fallbacks;
                        paths=decode_leaf_range(
                            channel_llrs,
                            std::move(paths),
                            frozen,
                            op.start,
                            op.size,
                            list_size,
                            stats);
                        continue;
                    }
                    node.candidates.resize(list_size);
                }
                stats.special_candidates+=node.stats.candidates_materialized;
                std::vector<DecoderPath> expanded;
                expanded.reserve(node.candidates.size());
                for(const Rate1Candidate& candidate:node.candidates)
                {
                    auto u=paths[candidate.parent_id].u;
                    const auto local_u=polar_transform(candidate.bits);
                    u.insert(u.end(),local_u.begin(),local_u.end());
                    expanded.push_back({candidate.path_metric,std::move(u)});
                }
                paths=prune_paths(std::move(expanded),list_size);
                ++stats.sorts;
                continue;
            }

            if(op.kind == DecodeOp::Kind::Rate1)
            {
                std::vector<ParentPath> parents;
                for(std::size_t i = 0; i < paths.size(); ++i)
                parents.push_back({i, paths[i].metric, node_alpha(channel_llrs, paths[i].u, op.start, op.size)});
                if(options.enable_rate1_tie_guard&&rate1_metric_tie_risk(parents))
                {
                    ++result_stats.tie_fallbacks;
                    paths=decode_leaf_range(channel_llrs,std::move(paths),frozen,op.start,op.size,list_size,stats);
                    continue;
                }
                auto choice = choose_rate1_ir(algorithm, parents, list_size);
                ++result_stats.engine_counts[choice.ir.engine];
                CandidateRankFn rank_fn = [&](const Rate1Candidate& candidate)
                {
                    auto local_u = polar_transform(candidate.bits);
                    std::vector<std::uint64_t> tie;
                    for(Bit b : paths[candidate.parent_id].u) tie.push_back(b);
                    for(Bit b : local_u) tie.push_back(b);
                    return GeneralRankKey{candidate.path_metric, std::move(tie)};
                };
                HiddenRankFn hidden_fn = [&](const ParentPath& parent, double metric)
                {
                    std::vector<std::uint64_t> tie;
                    for(Bit b : paths[parent.parent_id].u) tie.push_back(b);
                    tie.insert(tie.end(), op.size, 0);
                    return GeneralRankKey{metric, std::move(tie)};
                };
                auto integrated=execute_rate1_integrated(
                    choice.ir,
                    parents,
                    list_size,
                    rank_fn,
                    hidden_fn,
                    options.enable_rate1_tie_guard);
                stats.include_rate1(integrated.stats);
                if(!integrated.candidates || integrated.stats.fallback_used)
                {
                    if(integrated.stats.fallback_reason == "whole-frame tie at Top-L cutoff") ++result_stats.tie_fallbacks;
                    paths = decode_leaf_range(channel_llrs, std::move(paths), frozen, op.start, op.size, list_size, stats);
                    continue;
                }

                std::vector<DecoderPath> expanded;
                for(const auto& candidate : *integrated.candidates)
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

        if(options.enable_final_list_tie_guard&&has_close_path_metrics(paths))
        {
            ++result_stats.tie_fallbacks;
            auto canonical=scl_decode(channel_llrs,frozen,list_size);
            return {std::move(canonical.paths),std::move(result_stats)};
        }

        return {std::move(paths),std::move(result_stats)};
    }
}

// namespace polar_scl
