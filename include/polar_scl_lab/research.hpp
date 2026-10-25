#pragma once
#include "polar_scl_lab/linear.hpp"
#include "polar_scl/construction.hpp"
#include "polar_scl/search.hpp"
#include <array>
#include <map>
#include <random>

namespace polar_scl_lab {
using polar_scl::EvaluationConfig;using polar_scl::PolicyResult;using polar_scl::evaluate_policy;using polar_scl::pareto_policies;
using polar_scl::SearchSpace;using polar_scl::EvolutionConfig;using polar_scl::SearchResult;using polar_scl::evolutionary_search;
Rate1Policy random_policy(const SearchSpace&,std::mt19937_64&);
Rate1Policy mutate_policy(const Rate1Policy&,const SearchSpace&,std::mt19937_64&);
using polar_scl::polarization_weight_order;using polar_scl::pw_frozen_pattern;using polar_scl::place_information;
using polar_scl::extract_information;using polar_scl::append_crc;using polar_scl::crc_ok;using polar_scl::encode_payload;
using polar_scl::select_crc_path;using polar_scl::DecoderPath;using polar_scl::polar_transform;using polar_scl::scl_decode;

struct AlgorithmIR {
    std::string name,domain,ordering,generator,constraint,hidden_bound,termination,fallback;
    Rate1Policy policy{};
    std::string origin{"known"};std::vector<std::string> parents;std::string derivation_rule;
    bool exact_by_construction() const;
    std::string operator_signature() const;
    std::string structural_signature() const;
    std::string to_json() const;
};
struct PolarDecoderIR {AlgorithmIR rate1,spc,gpc2,gpc4;bool enable_spc=true,enable_gpc2=true,enable_gpc4=true;};
void validate_algorithm_ir(const AlgorithmIR&);
void validate_polar_decoder_ir(const PolarDecoderIR&);
std::vector<AlgorithmIR> known_algorithm_library(const Rate1Policy& = {0,1,ExpansionMode::MinimumBound,256});
struct DerivationTrace {std::string rule;std::vector<std::string> parents,inherited_operators,synthesized_operators,safety_obligations;};
struct DerivedAlgorithm {AlgorithmIR algorithm;DerivationTrace trace;};
struct OperatorPool {std::vector<std::string> domains,orderings,generators,constraints,hidden_bounds,terminations,fallbacks;};
struct SynthesisReport {
    OperatorPool operator_pool;std::size_t raw_compositions{},type_valid_compositions{},rejected_compositions{},known_operator_signatures{},novel_operator_signatures{};
    std::vector<DerivedAlgorithm> derived;
};
std::optional<std::string> semantic_hidden_bound_closure(const std::string&,const std::string&);
OperatorPool mine_operator_pool(const std::vector<AlgorithmIR>&);
SynthesisReport synthesize_algorithm_structures(std::vector<AlgorithmIR> = {}, std::vector<Rate1Policy> = {});
DerivedAlgorithm derive_constraint_transferred_frontier(const AlgorithmIR&,const AlgorithmIR&,std::string name = {},std::optional<Rate1Policy> = {});
std::vector<DerivedAlgorithm> derive_novel_algorithm_structures(std::vector<Rate1Policy> = {});
PolarDecoderIR default_polar_decoder_ir(const Rate1Policy& = {0,1,ExpansionMode::MinimumBound,256});
struct MaterializedResult {AlgorithmIR algorithm;std::optional<std::vector<Rate1Candidate>> candidates;SearchStats stats;};
MaterializedResult execute_algorithm_ir(const AlgorithmIR&,const std::vector<ParentPath>&,std::size_t,CandidateRankFn = {},HiddenRankFn = {});
MaterializedResult execute_algorithm_ir_integrated(const AlgorithmIR&,const std::vector<ParentPath>&,std::size_t,CandidateRankFn,HiddenRankFn);
struct IRStructureEvaluation {
    AlgorithmIR algorithm;std::size_t cases{},mismatches{};double mean_candidate_work{},p99_candidate_work{},fallback_rate{},mean_frontier_expansions{};
    bool is_exact()const{return mismatches==0;}
    std::string to_json()const;
};
struct IRStructureSearchResult {std::vector<IRStructureEvaluation> evaluated;IRStructureEvaluation selected,validation;std::vector<IRStructureEvaluation> selected_by_domain,validation_by_domain;};
IRStructureEvaluation evaluate_algorithm_structure(const AlgorithmIR&,std::size_t trials_per_case=25,std::uint64_t seed=20260829);
IRStructureSearchResult search_derived_algorithm_structures(std::size_t trials_per_case=25,std::uint64_t seed=20260829,std::uint64_t validation_seed=20270829);
std::string json_quote(const std::string&);
}
