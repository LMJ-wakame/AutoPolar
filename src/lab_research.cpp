#include "polar_scl_lab/research.hpp"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

namespace polar_scl_lab {
std::string linear_constraint_json(const LinearConstraintIR&h){std::ostringstream o;o<<"{\"name\":"<<json_quote(h.name)<<",\"num_variables\":"<<h.num_variables<<",\"declared_checks\":"<<h.parity_checks.size()<<",\"effective_rank\":"<<h.rank()<<",\"state_count\":"<<h.state_count()<<",\"parity_checks\":[";for(std::size_t i=0;i<h.parity_checks.size();++i){if(i)o<<',';o<<h.parity_checks[i];}o<<"],\"compiled_columns\":[";auto columns=h.columns();for(std::size_t i=0;i<columns.size();++i){if(i)o<<',';o<<columns[i];}o<<"]}";return o.str();}
std::string ConstraintProofMetadata::to_json()const{std::ostringstream o;o<<"{\"constraint_name\":"<<json_quote(constraint_name)<<",\"effective_rank\":"<<effective_rank<<",\"state_count\":"<<state_count<<",\"recurrence\":"<<json_quote(recurrence)<<",\"safety_obligations\":[";for(std::size_t i=0;i<safety_obligations.size();++i){if(i)o<<',';o<<json_quote(safety_obligations[i]);}o<<"]}";return o.str();}
std::string CompiledLinearAlgorithm::to_json()const{std::ostringstream o;const char*name=engine==LinearEngine::LazyExact?"lazy_exact":engine==LinearEngine::CertifiedFrontier?"certified_frontier":engine==LinearEngine::SyndromeAStar?"syndrome_astar":"sparse_trellis";o<<"{\"constraint\":"<<linear_constraint_json(constraint)<<",\"engine\":"<<json_quote(name)<<",\"state_budget\":";if(state_budget)o<<*state_budget;else o<<"null";o<<",\"proof_status\":\"obligations_emitted_not_formally_verified\",\"constraint_obligations\":"<<proof_metadata(constraint).to_json()<<",\"ordering\":\"metric,parent_id,mask\",\"fallback\":\"lazy_exact\"}";return o.str();}
Rate1Policy random_policy(const SearchSpace&s,std::mt19937_64&rng){s.validate();auto pick=[&](const auto&v){return v[std::uniform_int_distribution<std::size_t>(0,v.size()-1)(rng)];};return {pick(s.initial_active_bits),pick(s.growth_steps),pick(s.expansion_modes),pick(s.candidate_budgets)};}
Rate1Policy mutate_policy(const Rate1Policy&p,const SearchSpace&s,std::mt19937_64&rng) {
    s.validate(); std::vector<unsigned> fields;
    auto differs=[](const auto&v,auto current){return std::any_of(v.begin(),v.end(),[&](auto x){return x!=current;});};
    if(differs(s.initial_active_bits,p.initial_active_bits))fields.push_back(0);
    if(differs(s.growth_steps,p.growth_step))fields.push_back(1);
    if(differs(s.expansion_modes,p.expansion_mode))fields.push_back(2);
    if(differs(s.candidate_budgets,p.candidate_budget))fields.push_back(3);
    if(fields.empty())return p;
    auto choose=[&](const auto&v,auto current){auto choices=v;choices.erase(std::remove(choices.begin(),choices.end(),current),choices.end());return choices[std::uniform_int_distribution<std::size_t>(0,choices.size()-1)(rng)];};
    auto child=p;switch(fields[std::uniform_int_distribution<std::size_t>(0,fields.size()-1)(rng)]) {
        case 0:child.initial_active_bits=choose(s.initial_active_bits,p.initial_active_bits);break;
        case 1:child.growth_step=choose(s.growth_steps,p.growth_step);break;
        case 2:child.expansion_mode=choose(s.expansion_modes,p.expansion_mode);break;
        case 3:child.candidate_budget=choose(s.candidate_budgets,p.candidate_budget);break;
    }
    return child;
}
std::string json_quote(const std::string&s){std::ostringstream o;o<<'"';for(unsigned char c:s){switch(c){case '"':o<<"\\\"";break;case '\\':o<<"\\\\";break;case '\n':o<<"\\n";break;case '\r':o<<"\\r";break;case '\t':o<<"\\t";break;default:if(c<32)o<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<static_cast<unsigned>(c)<<std::dec;else o<<c;}}o<<'"';return o.str();}
bool AlgorithmIR::exact_by_construction()const{return (generator=="full_enumeration"&&termination=="enumerate_all")||(generator=="lazy_subset_stream"&&termination=="stream_top_l")||(generator=="bounded_frontier"&&termination=="hidden_bound_certificate"&&fallback=="lazy_exact");}
std::string AlgorithmIR::operator_signature()const{return domain+"|"+ordering+"|"+generator+"|"+constraint+"|"+hidden_bound+"|"+termination+"|"+fallback;}
std::string AlgorithmIR::structural_signature()const{return operator_signature()+"|"+polar_scl::to_string(policy.expansion_mode);}
std::string AlgorithmIR::to_json()const {std::ostringstream o;o<<"{\"name\":"<<json_quote(name)<<",\"domain\":"<<json_quote(domain)<<",\"operators\":{\"ordering\":"<<json_quote(ordering)<<",\"generator\":"<<json_quote(generator)<<",\"constraint\":"<<json_quote(constraint)<<",\"hidden_bound\":"<<json_quote(hidden_bound)<<",\"termination\":"<<json_quote(termination)<<",\"fallback\":"<<json_quote(fallback)<<"},\"policy\":{\"initial_active_bits\":"<<policy.initial_active_bits<<",\"growth_step\":"<<policy.growth_step<<",\"expansion_mode\":"<<json_quote(polar_scl::to_string(policy.expansion_mode))<<",\"candidate_budget\":"<<policy.candidate_budget<<"},\"origin\":"<<json_quote(origin)<<",\"parents\":[";for(std::size_t i=0;i<parents.size();++i){if(i)o<<',';o<<json_quote(parents[i]);}o<<"],\"derivation_rule\":"<<json_quote(derivation_rule)<<",\"exact_by_construction\":"<<(exact_by_construction()?"true":"false")<<'}';return o.str();}
namespace {
const std::map<std::string,std::string> constraints{{"rate1","unconstrained"},{"spc","even_parity"},{"gpc2","two_even_parities"},{"gpc4","four_even_parities"}};
std::string bound_for(const std::string&c){if(c=="unconstrained")return "next_hidden_reliability";if(c=="even_parity")return "parity_completion";if(c=="two_even_parities"||c=="four_even_parities")return "syndrome_completion";throw std::invalid_argument("unknown constraint");}
DerivationTrace trace_for(const AlgorithmIR&a){return {a.derivation_rule,a.parents,{"ordering,generator,termination,fallback inherited from frontier parent","constraint inherited from exact parent"},
    {a.hidden_bound=="parity_completion"?"min(one-hidden plus active parity repair,two-hidden plus active parity repair)":"min-plus XOR convolution of active and nonempty-hidden syndrome tables"},
    {"visible candidates satisfy constraint","hidden bound covers all unseen feasible candidates","strict certificate","budget invokes exact fallback"}};}
bool same(const std::vector<Rate1Candidate>&a,const std::vector<Rate1Candidate>&b){if(a.size()!=b.size())return false;for(std::size_t i=0;i<a.size();++i)if(polar_scl::local_rank_key(a[i])!=polar_scl::local_rank_key(b[i]))return false;return true;}
}
void validate_algorithm_ir(const AlgorithmIR&a){validate_policy(a.policy);if(a.name.find_first_not_of(" \t\r\n")==std::string::npos)throw std::invalid_argument("empty IR name");
    auto expected=constraints.find(a.domain);if(expected==constraints.end()||a.constraint!=expected->second)throw std::invalid_argument("IR domain/constraint mismatch");
    std::array<std::string,4> want;
    if(a.generator=="full_enumeration")want={"natural_mask","none","enumerate_all","none"};
    else if(a.generator=="lazy_subset_stream")want={"reliability_ascending","none","stream_top_l","none"};
    else if(a.generator=="bounded_frontier")want={"reliability_ascending",bound_for(a.constraint),"hidden_bound_certificate","lazy_exact"};
    else throw std::invalid_argument("unsupported IR generator");
    if(want!=std::array<std::string,4>{a.ordering,a.hidden_bound,a.termination,a.fallback})throw std::invalid_argument("incompatible IR operators");
    if(a.origin=="known"&&(!a.parents.empty()||!a.derivation_rule.empty()))throw std::invalid_argument("known IR cannot carry derived provenance");
    if(a.origin=="derived"&&(a.parents.size()<2||a.derivation_rule.empty()))throw std::invalid_argument("derived IR requires parents and rule");
    if((a.origin!="known"&&a.origin!="derived")||!a.exact_by_construction())throw std::invalid_argument("invalid origin or unsafe composition");}
void validate_polar_decoder_ir(const PolarDecoderIR&d){const AlgorithmIR*slots[]={&d.rate1,&d.spc,&d.gpc2,&d.gpc4};const char*names[]={"rate1","spc","gpc2","gpc4"};for(std::size_t i=0;i<4;++i){validate_algorithm_ir(*slots[i]);if(slots[i]->domain!=names[i])throw std::invalid_argument("wrong decoder IR slot");}}
std::vector<AlgorithmIR> known_algorithm_library(const Rate1Policy&policy){std::vector<AlgorithmIR> lib;
    for(const std::string domain:{"rate1","spc","gpc2","gpc4"})for(const std::string gen:{"full_enumeration","lazy_subset_stream","bounded_frontier"}){
        if(gen=="bounded_frontier"&&domain!="rate1")continue;
        AlgorithmIR a;a.domain=domain;a.constraint=constraints.at(domain);a.generator=gen;
        const bool exhaustive=gen=="full_enumeration",front=gen=="bounded_frontier";
        a.name="known_"+domain+(exhaustive?"_exhaustive":front?"_certified_frontier":"_lazy_exact");a.ordering=exhaustive?"natural_mask":"reliability_ascending";
        a.hidden_bound=front?bound_for(a.constraint):"none";a.termination=exhaustive?"enumerate_all":front?"hidden_bound_certificate":"stream_top_l";a.fallback=front?"lazy_exact":"none";if(front)a.policy=policy;validate_algorithm_ir(a);lib.push_back(a);}
    return lib;}
std::optional<std::string> semantic_hidden_bound_closure(const std::string&generator,const std::string&constraint){if(generator!="bounded_frontier")return "none";try{return bound_for(constraint);}catch(const std::invalid_argument&){return {};}}
OperatorPool mine_operator_pool(const std::vector<AlgorithmIR>&lib){if(lib.empty())throw std::invalid_argument("empty IR library");std::array<std::set<std::string>,7> sets;
    for(const auto&a:lib){validate_algorithm_ir(a);std::array<std::string,7> values{a.domain,a.ordering,a.generator,a.constraint,a.hidden_bound,a.termination,a.fallback};for(std::size_t i=0;i<7;++i)sets[i].insert(values[i]);}
    for(const auto&g:sets[2])for(const auto&c:sets[3])if(auto b=semantic_hidden_bound_closure(g,c))sets[4].insert(*b);
    OperatorPool p;std::vector<std::string>*v[]={&p.domains,&p.orderings,&p.generators,&p.constraints,&p.hidden_bounds,&p.terminations,&p.fallbacks};for(std::size_t i=0;i<7;++i)v[i]->assign(sets[i].begin(),sets[i].end());return p;}
SynthesisReport synthesize_algorithm_structures(std::vector<AlgorithmIR>lib,std::vector<Rate1Policy>policies){if(lib.empty())lib=known_algorithm_library();if(policies.empty())policies={{0,1,ExpansionMode::MinimumBound,64},{0,1,ExpansionMode::MinimumBound,256},{1,1,ExpansionMode::MinimumBound,256},{0,2,ExpansionMode::AllCompetitive,256}};
    SynthesisReport r;r.operator_pool=mine_operator_pool(lib);const auto&p=r.operator_pool;std::set<std::string>known;for(const auto&a:lib)known.insert(a.operator_signature());r.known_operator_signatures=known.size();std::map<std::string,AlgorithmIR>valid;
    for(const auto&d:p.domains)for(const auto&o:p.orderings)for(const auto&g:p.generators)for(const auto&c:p.constraints)for(const auto&b:p.hidden_bounds)for(const auto&t:p.terminations)for(const auto&f:p.fallbacks){++r.raw_compositions;AlgorithmIR a;a.name="grammar_candidate";a.domain=d;a.ordering=o;a.generator=g;a.constraint=c;a.hidden_bound=b;a.termination=t;a.fallback=f;
        try{validate_algorithm_ir(a);valid[a.operator_signature()]=a;}catch(const std::invalid_argument&){++r.rejected_compositions;}}
    r.type_valid_compositions=valid.size();std::size_t structure=0;
    for(const auto&[signature,base]:valid){if(known.contains(signature))continue;++r.novel_operator_signatures;const AlgorithmIR*front=nullptr,*exact=nullptr;
        for(const auto&a:lib){if(!front&&a.generator==base.generator&&a.ordering==base.ordering&&a.termination==base.termination&&a.fallback==base.fallback)front=&a;
            if(a.domain==base.domain&&a.constraint==base.constraint&&(!exact||a.generator=="lazy_subset_stream"))exact=&a;}
        if(!front||!exact)continue;
        for(std::size_t i=0;i<policies.size();++i){auto a=base;a.name="derived_"+a.domain+"_"+a.generator+"_"+std::to_string(structure)+"_"+std::to_string(i);a.policy=policies[i];a.origin="derived";a.parents={front->name,exact->name};a.derivation_rule="operator_pool_enumeration_with_semantic_closure";validate_algorithm_ir(a);r.derived.push_back({a,trace_for(a)});}++structure;}
    return r;}
DerivedAlgorithm derive_constraint_transferred_frontier(const AlgorithmIR&front,const AlgorithmIR&exact,std::string name,std::optional<Rate1Policy>policy){validate_algorithm_ir(front);validate_algorithm_ir(exact);if(front.domain!="rate1"||front.generator!="bounded_frontier"||exact.domain=="rate1")throw std::invalid_argument("invalid constraint-transfer parents");
    auto a=front;a.name=name.empty()?"derived_"+exact.domain+"_certified_frontier":name;a.domain=exact.domain;a.constraint=exact.constraint;a.hidden_bound=bound_for(a.constraint);a.origin="derived";a.parents={front.name,exact.name};a.derivation_rule="constraint_transfer_with_bound_closure";a.policy=policy.value_or(front.policy);validate_algorithm_ir(a);return {a,trace_for(a)};}
std::vector<DerivedAlgorithm> derive_novel_algorithm_structures(std::vector<Rate1Policy>p){return synthesize_algorithm_structures({},p).derived;}
PolarDecoderIR default_polar_decoder_ir(const Rate1Policy&p){auto lib=known_algorithm_library(p);PolarDecoderIR d;for(const auto&a:lib)if(a.name=="known_rate1_certified_frontier")d.rate1=a;auto synthesis=synthesize_algorithm_structures(lib,{p});for(const auto&x:synthesis.derived){if(x.algorithm.domain=="spc")d.spc=x.algorithm;else if(x.algorithm.domain=="gpc2")d.gpc2=x.algorithm;else if(x.algorithm.domain=="gpc4")d.gpc4=x.algorithm;}validate_polar_decoder_ir(d);return d;}
MaterializedResult execute_algorithm_ir(const AlgorithmIR&a,const std::vector<ParentPath>&p,std::size_t L,CandidateRankFn rank,HiddenRankFn hidden){validate_algorithm_ir(a);Rate1SearchResult r;
    if(a.domain=="rate1"){if(a.generator=="full_enumeration")r=exhaustive_rate1_top_l(p,L);else if(a.generator=="lazy_subset_stream")r=lazy_exact_rate1_top_l(p,L);else r=adaptive_rate1_top_l(p,L,a.policy,rank,hidden);}
    else if(a.domain=="spc"){if(a.generator=="full_enumeration")r=exhaustive_spc_top_l(p,L);else if(a.generator=="lazy_subset_stream")r=lazy_exact_spc_top_l(p,L);else r=adaptive_spc_top_l(p,L,a.policy,rank,hidden);}
    else {auto checks=a.domain=="gpc2"?2u:4u;if(a.generator=="full_enumeration")r=exhaustive_gpc_top_l(p,L,checks);else if(a.generator=="lazy_subset_stream")r=lazy_exact_gpc_top_l(p,L,checks);else r=adaptive_gpc_top_l(p,L,checks,a.policy,rank,hidden);}
    return {a,r.candidates,r.stats};}
MaterializedResult execute_algorithm_ir_integrated(const AlgorithmIR&a,const std::vector<ParentPath>&p,std::size_t L,CandidateRankFn rank,HiddenRankFn hidden){validate_algorithm_ir(a);if(a.generator!="bounded_frontier"){SearchStats s;s.algorithm="ir_"+a.name+"_canonical_fallback";s.fallback_used=true;s.fallback_reason="IR engine lacks whole-prefix tie integration";return {a,{},s};}return execute_algorithm_ir(a,p,L,rank,hidden);}
std::string IRStructureEvaluation::to_json()const{std::ostringstream o;o<<std::setprecision(17)<<"{\"algorithm\":"<<algorithm.to_json()<<",\"cases\":"<<cases<<",\"mismatches\":"<<mismatches<<",\"mean_candidate_work\":"<<mean_candidate_work<<",\"p99_candidate_work\":"<<p99_candidate_work<<",\"fallback_rate\":"<<fallback_rate<<",\"mean_frontier_expansions\":"<<mean_frontier_expansions<<'}';return o.str();}
IRStructureEvaluation evaluate_algorithm_structure(const AlgorithmIR&a,std::size_t trials,std::uint64_t seed){if(trials==0)throw std::invalid_argument("trials must be positive");validate_algorithm_ir(a);std::mt19937_64 rng(seed);std::uniform_int_distribution<int>base(0,8),llr(-16,16);IRStructureEvaluation r;r.algorithm=a;std::vector<std::size_t>work;double expansions=0,fallbacks=0;
    auto oracle=a;oracle.generator="lazy_subset_stream";oracle.ordering="reliability_ascending";oracle.hidden_bound="none";oracle.termination="stream_top_l";oracle.fallback="none";
    for(std::size_t n:{8u,16u,32u})for(std::size_t L:{2u,4u,8u})for(std::size_t k=0;k<trials;++k){std::vector<ParentPath>parents;for(std::size_t id=0;id<L;++id){ParentPath p{id,base(rng)/4.0,{}};for(std::size_t i=0;i<n;++i)p.llrs.push_back(llr(rng)/4.0);parents.push_back(p);}
        auto ex=execute_algorithm_ir(oracle,parents,L),actual=execute_algorithm_ir(a,parents,L);++r.cases;r.mismatches+=!same(*ex.candidates,*actual.candidates);work.push_back(actual.stats.candidates_materialized);expansions+=static_cast<double>(actual.stats.frontier_expansions);fallbacks+=actual.stats.fallback_used?1:0;}
    r.mean_candidate_work=std::accumulate(work.begin(),work.end(),0.0)/static_cast<double>(r.cases);std::sort(work.begin(),work.end());r.p99_candidate_work=static_cast<double>(work[std::min(work.size()-1,static_cast<std::size_t>(.99*static_cast<double>(work.size())))]);r.fallback_rate=fallbacks/static_cast<double>(r.cases);r.mean_frontier_expansions=expansions/static_cast<double>(r.cases);return r;}
IRStructureSearchResult search_derived_algorithm_structures(std::size_t trials,std::uint64_t seed,std::uint64_t validation_seed){if(seed==validation_seed)throw std::invalid_argument("training and validation seeds must differ");IRStructureSearchResult r;for(const auto&x:derive_novel_algorithm_structures())r.evaluated.push_back(evaluate_algorithm_structure(x.algorithm,trials,seed));
    std::stable_sort(r.evaluated.begin(),r.evaluated.end(),[](const auto&a,const auto&b){return std::tuple{!a.is_exact(),a.mean_candidate_work,a.p99_candidate_work,a.fallback_rate}<std::tuple{!b.is_exact(),b.mean_candidate_work,b.p99_candidate_work,b.fallback_rate};});
    if(r.evaluated.empty()||!r.evaluated.front().is_exact()){throw std::runtime_error("no exact IR found");}
    r.selected=r.evaluated.front();
    for(const std::string domain:{"gpc2","gpc4","spc"}){auto found=std::find_if(r.evaluated.begin(),r.evaluated.end(),[&](const auto&e){return e.algorithm.domain==domain;});if(found==r.evaluated.end())throw std::runtime_error("missing derived domain");r.selected_by_domain.push_back(*found);auto check=evaluate_algorithm_structure(found->algorithm,trials,validation_seed);if(!check.is_exact())throw std::runtime_error("independent validation failed");r.validation_by_domain.push_back(check);if(domain==r.selected.algorithm.domain)r.validation=check;}return r;}
}
