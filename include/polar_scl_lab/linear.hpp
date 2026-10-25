#pragma once
#include "polar_scl/linear_trellis.hpp"
#include <optional>

namespace polar_scl_lab {
using polar_scl::Bit;
using polar_scl::ParentPath;
using polar_scl::Rate1Candidate;
using polar_scl::SearchStats;
using polar_scl::Rate1SearchResult;
using polar_scl::Rate1Policy;
using polar_scl::ExpansionMode;
using polar_scl::CandidateRankFn;
using polar_scl::HiddenRankFn;
using polar_scl::LinearConstraintIR;
using polar_scl::TrellisIR;
using polar_scl::TrellisSearchReport;
using polar_scl::search_trellis_ir;
using polar_scl::linear_trellis_top_l;
using polar_scl::exhaustive_rate1_top_l;
using polar_scl::lazy_exact_rate1_top_l;
using polar_scl::adaptive_rate1_top_l;
using polar_scl::validate_policy;

struct ConstraintProofMetadata {
    std::string constraint_name, recurrence;
    std::size_t effective_rank{}, state_count{};
    std::vector<std::string> safety_obligations;
    std::string to_json() const;
};
std::string linear_constraint_json(const LinearConstraintIR&);
ConstraintProofMetadata proof_metadata(const LinearConstraintIR&);
LinearConstraintIR gpc_constraint(std::size_t n, std::size_t checks);
Rate1Candidate make_candidate(const ParentPath&, std::uint64_t);
double metric_sum(const std::vector<double>&);
Rate1SearchResult exhaustive_linear_top_l(const std::vector<ParentPath>&, std::size_t, const LinearConstraintIR&);
Rate1SearchResult lazy_exact_linear_top_l(const std::vector<ParentPath>&, std::size_t, const LinearConstraintIR&);
Rate1SearchResult adaptive_linear_top_l(const std::vector<ParentPath>&, std::size_t, const LinearConstraintIR&,
    const Rate1Policy& = {}, CandidateRankFn = {}, HiddenRankFn = {});
Rate1SearchResult astar_linear_top_l(const std::vector<ParentPath>&, std::size_t, const LinearConstraintIR&,
    std::optional<std::size_t> state_budget = {});
Rate1SearchResult exhaustive_spc_top_l(const std::vector<ParentPath>&, std::size_t);
Rate1SearchResult lazy_exact_spc_top_l(const std::vector<ParentPath>&, std::size_t);
Rate1SearchResult adaptive_spc_top_l(const std::vector<ParentPath>&, std::size_t, const Rate1Policy& = {}, CandidateRankFn = {}, HiddenRankFn = {});
Rate1SearchResult exhaustive_gpc_top_l(const std::vector<ParentPath>&, std::size_t, std::size_t);
Rate1SearchResult lazy_exact_gpc_top_l(const std::vector<ParentPath>&, std::size_t, std::size_t);
Rate1SearchResult adaptive_gpc_top_l(const std::vector<ParentPath>&, std::size_t, std::size_t,
    const Rate1Policy& = {}, CandidateRankFn = {}, HiddenRankFn = {});
inline Rate1SearchResult exhaustive_gpc2_top_l(const std::vector<ParentPath>&p,std::size_t L){return exhaustive_gpc_top_l(p,L,2);}
inline Rate1SearchResult exhaustive_gpc4_top_l(const std::vector<ParentPath>&p,std::size_t L){return exhaustive_gpc_top_l(p,L,4);}
inline Rate1SearchResult lazy_exact_gpc2_top_l(const std::vector<ParentPath>&p,std::size_t L){return lazy_exact_gpc_top_l(p,L,2);}
inline Rate1SearchResult lazy_exact_gpc4_top_l(const std::vector<ParentPath>&p,std::size_t L){return lazy_exact_gpc_top_l(p,L,4);}
inline Rate1SearchResult adaptive_gpc2_top_l(const std::vector<ParentPath>&p,std::size_t L,const Rate1Policy&policy={},CandidateRankFn rank={},HiddenRankFn hidden={}){return adaptive_gpc_top_l(p,L,2,policy,rank,hidden);}
inline Rate1SearchResult adaptive_gpc4_top_l(const std::vector<ParentPath>&p,std::size_t L,const Rate1Policy&policy={},CandidateRankFn rank={},HiddenRankFn hidden={}){return adaptive_gpc_top_l(p,L,4,policy,rank,hidden);}
enum class LinearEngine { LazyExact, CertifiedFrontier, SyndromeAStar, SparseTrellis };
struct CompiledLinearAlgorithm {
    LinearConstraintIR constraint;
    LinearEngine engine{LinearEngine::SyndromeAStar};
    Rate1Policy policy{};
    std::optional<std::size_t> state_budget;
    CompiledLinearAlgorithm(LinearConstraintIR,LinearEngine = LinearEngine::SyndromeAStar,
        std::optional<Rate1Policy> = {},std::optional<std::size_t> = {});
    Rate1SearchResult execute(const std::vector<ParentPath>&, std::size_t) const;
    std::string to_json() const;
};
CompiledLinearAlgorithm compile_linear_algorithm_ir(const LinearConstraintIR&, LinearEngine = LinearEngine::SyndromeAStar,
    std::optional<Rate1Policy> = {}, std::optional<std::size_t> = {});
}
