#include "test_util.hpp"
#include "polar_scl_lab/polar_scl_lab.hpp"
#include <random>
#include <set>

namespace lab=polar_scl_lab;
int main(){return run_test_main("full_cpp_migration",[]{
    std::size_t node_cases=0,frame_cases=0,groups=0;
    auto group=[&](const char*name){++groups;std::cout<<"[section] "<<name<<'\n';};
    std::mt19937_64 rng(701);std::uniform_int_distribution<int> w(-5,5),row(1,255);
    auto parents=[&](std::size_t n,std::size_t count){std::vector<lab::ParentPath>out;for(std::size_t p=0;p<count;++p){lab::ParentPath x{p,static_cast<double>(p),{}};for(std::size_t i=0;i<n;++i)x.llrs.push_back(w(rng));out.push_back(x);}return out;};
    group("rank compression and parity equivalence");
    lab::LinearConstraintIR redundant{"redundant",4,{3,12,15}};CHECK_EQ(redundant.rank(),2u);CHECK_EQ(redundant.state_count(),4u);
    for(std::uint64_t a=0;a<16;++a)for(std::uint64_t b=0;b<16;++b){bool equal=true;for(auto h:redundant.parity_checks)equal=equal&&((std::popcount(a&h)&1)==(std::popcount(b&h)&1));CHECK_EQ(redundant.syndrome(a)==redundant.syndrome(b),equal);}
    CHECK_THROW((lab::LinearConstraintIR{"bad",4,{16}}).validate());
    CHECK_THROW((lab::LinearConstraintIR{"bad",0,{1}}).validate());
    group("arbitrary-H exhaustive/lazy/frontier/A-star/compiler/trellis");
    for(int matrix=0;matrix<12;++matrix){lab::LinearConstraintIR h{"random",8,{static_cast<std::uint64_t>(row(rng)),static_cast<std::uint64_t>(row(rng)),static_cast<std::uint64_t>(row(rng))}};
        auto trellis=lab::search_trellis_ir(h,32,73);
        for(std::size_t L:{1u,2u,4u,8u})for(int trial=0;trial<5;++trial){auto p=parents(8,2);auto ex=lab::exhaustive_linear_top_l(p,L,h);
            check_candidates_equal(ex.candidates,lab::lazy_exact_linear_top_l(p,L,h).candidates);
            check_candidates_equal(ex.candidates,lab::adaptive_linear_top_l(p,L,h,{0,1,lab::ExpansionMode::MinimumBound,512}).candidates);
            check_candidates_equal(ex.candidates,lab::adaptive_linear_top_l(p,L,h,{0,2,lab::ExpansionMode::AllCompetitive,512}).candidates);
            check_candidates_equal(ex.candidates,lab::astar_linear_top_l(p,L,h).candidates);
            check_candidates_equal(ex.candidates,lab::linear_trellis_top_l(p,L,trellis).candidates);
            for(auto engine:{lab::LinearEngine::LazyExact,lab::LinearEngine::CertifiedFrontier,lab::LinearEngine::SyndromeAStar})check_candidates_equal(ex.candidates,lab::compile_linear_algorithm_ir(h,engine).execute(p,L).candidates);
            ++node_cases;}}
    group("exact budgets and counters");
    auto p=parents(8,4);auto h=lab::gpc_constraint(8,2);auto ex=lab::exhaustive_linear_top_l(p,8,h);
    auto astar=lab::astar_linear_top_l(p,8,h,2);CHECK_TRUE(astar.stats.fallback_used);check_candidates_equal(ex.candidates,astar.candidates);
    auto front=lab::adaptive_linear_top_l(p,8,h,{0,1,lab::ExpansionMode::MinimumBound,1});CHECK_TRUE(front.stats.fallback_used);check_candidates_equal(ex.candidates,front.candidates);
    auto counted=lab::adaptive_linear_top_l(p,8,h,{0,1,lab::ExpansionMode::MinimumBound,4096});CHECK_EQ(counted.stats.dp_state_relaxations,4u*41u*4u);CHECK_TRUE(counted.stats.primitive_work()>0);
    CHECK_THROW(lab::compile_linear_algorithm_ir(h,lab::LinearEngine::LazyExact,{},2));
    CHECK_THROW(lab::compile_linear_algorithm_ir(h,lab::LinearEngine::LazyExact,lab::Rate1Policy{}));
    auto wide=parents(64,2);auto wide_h=lab::gpc_constraint(64,2);
    auto wide_ref=lab::lazy_exact_linear_top_l(wide,4,wide_h);
    check_candidates_equal(wide_ref.candidates,lab::adaptive_linear_top_l(wide,4,wide_h).candidates);
    check_candidates_equal(wide_ref.candidates,lab::astar_linear_top_l(wide,4,wide_h,10000).candidates);
    group("SPC and generic GPC2/4/8");
    for(std::size_t checks:{1u,2u,4u,8u})for(std::size_t L:{1u,2u,4u,8u})for(int trial=0;trial<15;++trial){auto in=parents(checks==8?16:8,L);
        auto ref=checks==1?lab::exhaustive_spc_top_l(in,L):lab::exhaustive_gpc_top_l(in,L,checks);
        auto lazy=checks==1?lab::lazy_exact_spc_top_l(in,L):lab::lazy_exact_gpc_top_l(in,L,checks);
        auto adaptive=checks==1?lab::adaptive_spc_top_l(in,L):lab::adaptive_gpc_top_l(in,L,checks);
        check_candidates_equal(ref.candidates,lazy.candidates);check_candidates_equal(ref.candidates,adaptive.candidates);++node_cases;}
    group("typed grammar and derivation provenance");
    auto grammar=lab::synthesize_algorithm_structures();CHECK_EQ(grammar.raw_compositions,2304u);CHECK_EQ(grammar.type_valid_compositions,12u);CHECK_EQ(grammar.rejected_compositions,2292u);CHECK_EQ(grammar.novel_operator_signatures,3u);CHECK_EQ(grammar.derived.size(),12u);
    auto known=lab::known_algorithm_library();CHECK_EQ(known.size(),9u);auto derived=lab::derive_novel_algorithm_structures({{0,1,lab::ExpansionMode::MinimumBound,256}});CHECK_EQ(derived.size(),3u);
    auto wrong=derived.front().algorithm;wrong.fallback="none";CHECK_THROW(lab::validate_algorithm_ir(wrong));
    for(const auto&a:known){lab::validate_algorithm_ir(a);CHECK_TRUE(a.exact_by_construction());}
    for(const auto&d:derived){CHECK_EQ(d.trace.parents.size(),2u);CHECK_TRUE(!d.trace.synthesized_operators.empty());}
    auto all=known;for(const auto&d:derived)all.push_back(d.algorithm);
    for(const auto&a:all)for(int trial=0;trial<20;++trial){auto in=parents(8,4);auto oracle=a;oracle.generator="full_enumeration";oracle.ordering="natural_mask";oracle.hidden_bound="none";oracle.termination="enumerate_all";oracle.fallback="none";
        auto ref=lab::execute_algorithm_ir(oracle,in,4),actual=lab::execute_algorithm_ir(a,in,4);check_candidates_equal(*ref.candidates,*actual.candidates);++node_cases;}
    auto lazy_ir=known[1];auto material=lab::execute_algorithm_ir_integrated(lazy_ir,p,4,polar_scl::general_local_rank_key,[](const lab::ParentPath&x,double m){return polar_scl::GeneralRankKey{m,{x.parent_id,0}};});CHECK_TRUE(!material.candidates.has_value());CHECK_TRUE(material.stats.fallback_used);
    group("schedule and risk-aware fusion/decomposition");
    std::vector<bool>flags{true,true,false,false,true,false,false,false};auto ops=lab::rate1_schedule(flags);CHECK_EQ(ops.size(),2u);CHECK_EQ(ops[0].kind,"gpc2");CHECK_EQ(ops[1].kind,"spc");
    std::vector<bool>g4flags(16,false);std::fill_n(g4flags.begin(),4,true);
    lab::ScheduleCostModel expensive;expensive.profiles={{"gpc4",16,8,1000,1000,0}};auto split=lab::risk_aware_schedule(g4flags,8,expensive);CHECK_TRUE(split.operations.size()>1);
    auto cheap=expensive;cheap.profiles[0].mean_work=cheap.profiles[0].p99_work=.1;auto fused=lab::risk_aware_schedule(g4flags,8,cheap);CHECK_EQ(fused.operations.size(),1u);CHECK_EQ(fused.operations[0].kind,"gpc4");
    group("whole-frame continuous and quantized canonical lists");
    std::normal_distribution<double>normal(0,3);
    for(std::size_t n:{4u,8u,16u,32u})for(std::size_t L:{1u,2u,4u,8u})for(int trial=0;trial<30;++trial){std::vector<double>llrs(n);std::vector<bool>frozen(n);for(std::size_t i=0;i<n;++i){llrs[i]=trial%2==0?static_cast<double>(w(rng)):normal(rng);frozen[i]=(rng()%2)==0;}
        auto expected=lab::scl_decode(llrs,frozen,L);for(const lab::Rate1Policy policy:std::vector<lab::Rate1Policy>{{0,1,lab::ExpansionMode::MinimumBound,64},{0,2,lab::ExpansionMode::AllCompetitive,256},{0,1,lab::ExpansionMode::MinimumBound,1}}){auto actual=lab::rate1_scl_decode(llrs,frozen,L,policy);check_paths_equal(actual.paths,expected.paths);}
        ++frame_cases;}
    for(int trial=0;trial<40;++trial){std::vector<double>llrs(16);for(auto&v:llrs)v=normal(rng);auto ref=lab::scl_decode(llrs,g4flags,8);check_paths_equal(lab::rate1_scl_decode(llrs,g4flags,8,{},true,true,true,nullptr,&expensive).paths,ref.paths);check_paths_equal(lab::rate1_scl_decode(llrs,g4flags,8,{},true,true,true,nullptr,&cheap).paths,ref.paths);++frame_cases;}
    group("PW CRC and encoding");auto frozen=lab::pw_frozen_pattern(64,40);std::vector<lab::Bit>payload(32);for(auto&b:payload)b=static_cast<lab::Bit>(rng()%2);auto [u,x]=lab::encode_payload(payload,frozen);CHECK_EQ(lab::polar_transform(x),u);CHECK_TRUE(lab::crc_ok(lab::extract_information(u,frozen)));auto crc=lab::append_crc(payload);crc[0]^=1;CHECK_TRUE(!lab::crc_ok(crc));
    group("policy evolution and independent IR validation");lab::EvaluationConfig cfg;cfg.node_sizes={8};cfg.list_sizes={2};cfg.trials_per_case=2;auto evaluated=lab::evaluate_policy({},cfg);CHECK_TRUE(evaluated.is_exact());lab::EvolutionConfig evolution;evolution.generations=1;evolution.population_size=2;evolution.offspring_per_generation=2;auto evolution_result=lab::evolutionary_search({},evolution,cfg);CHECK_TRUE(!evolution_result.pareto_front.empty());
    lab::Rate1Policy parent_policy{0,1,lab::ExpansionMode::MinimumBound,32};auto child_policy=lab::mutate_policy(parent_policy,{},rng);
    CHECK_EQ((parent_policy.initial_active_bits!=child_policy.initial_active_bits)+(parent_policy.growth_step!=child_policy.growth_step)+(parent_policy.expansion_mode!=child_policy.expansion_mode)+(parent_policy.candidate_budget!=child_policy.candidate_budget),1);
    auto selected=lab::search_derived_algorithm_structures(1);CHECK_EQ(selected.selected_by_domain.size(),3u);for(const auto&v:selected.validation_by_domain)CHECK_TRUE(v.is_exact());CHECK_THROW(lab::search_derived_algorithm_structures(1,1,1));
    group("fractional parent metrics and equal-cost A-star terminals");
    std::exponential_distribution<double>base_metric(1);
    for(std::size_t L:{1u,2u,4u,8u})for(int sample=0;sample<50;++sample){auto inputs=parents(8,L);for(auto&x:inputs)x.path_metric=base_metric(rng);lab::LinearConstraintIR sh{"spc",8,{255}};auto oracle=lab::exhaustive_linear_top_l(inputs,L,sh);check_candidates_equal(oracle.candidates,lab::astar_linear_top_l(inputs,L,sh).candidates);check_candidates_equal(oracle.candidates,lab::adaptive_spc_top_l(inputs,L).candidates);++node_cases;}
    std::cout<<"migration_groups="<<groups<<" node_oracle_cases="<<node_cases<<" frame_cases="<<frame_cases<<'\n';
});}
