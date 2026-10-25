#pragma once
#include "polar_scl_lab/research.hpp"

namespace polar_scl_lab {
struct DecodeOp {std::string kind;std::size_t start{},size{};};
struct NodeCostEstimate {std::string kind;std::size_t size{},list_size{};double mean_work{},p99_work{},fallback_rate{};void validate()const;};
struct ScheduleCostModel {
    std::vector<NodeCostEstimate> profiles;
    double tail_weight=.25,fallback_penalty=100,frozen_leaf_cost=1,information_leaf_cost=2;
    void validate()const;
    std::optional<double> score(const std::string&,std::size_t,std::size_t)const;
};
struct ScheduleDecision {std::vector<DecodeOp> operations;double estimated_cost{};std::vector<std::string> rejected_constituents;};
struct DecoderStats {
    std::size_t leaf_steps{},path_splits{},sorts{},rate0_nodes{};
    std::size_t rate1_nodes{},rate1_candidates{},rate1_frontier_expansions{},rate1_fallbacks{};
    std::size_t spc_nodes{},spc_candidates{},spc_frontier_expansions{},spc_fallbacks{},spc_parity_rejections{};
    std::size_t gpc2_nodes{},gpc2_candidates{},gpc2_frontier_expansions{},gpc2_fallbacks{},gpc2_syndrome_rejections{};
    std::size_t gpc4_nodes{},gpc4_candidates{},gpc4_frontier_expansions{},gpc4_fallbacks{},gpc4_syndrome_rejections{};
    void include(const std::string&,const SearchStats&);
};
struct DecodeResult {std::vector<DecoderPath> paths;DecoderStats stats;};
std::vector<DecodeOp> rate1_schedule(const std::vector<bool>&,bool enable_spc=true,bool enable_gpc2=true,bool enable_gpc4=true);
ScheduleDecision risk_aware_schedule(const std::vector<bool>&,std::size_t,const ScheduleCostModel&,
    bool enable_spc=true,bool enable_gpc2=true,bool enable_gpc4=true);
DecodeResult rate1_scl_decode(const std::vector<double>&,const std::vector<bool>&,std::size_t,
    const Rate1Policy& = {0,1,ExpansionMode::MinimumBound,256},
    bool enable_spc=true,bool enable_gpc2=true,bool enable_gpc4=true,
    const PolarDecoderIR* decoder_ir=nullptr,const ScheduleCostModel* cost_model=nullptr);
}
