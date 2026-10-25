#include "polar_scl_lab/decoder.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <set>
#include <stdexcept>

namespace polar_scl_lab {
namespace {
void validate_flags(const std::vector<bool>&f){if(f.empty()||(f.size()&(f.size()-1)))throw std::invalid_argument("frozen pattern must have power-of-two length");}
std::string kind_at(const std::vector<bool>&f,std::size_t start,std::size_t n,bool spc,bool g2,bool g4){
    bool all=true,any=false;for(std::size_t i=start;i<start+n;++i){all=all&&f[i];any=any||f[i];}
    if(n>1&&all){return "rate0";}
    if(n>1&&!any){return "rate1";}
    auto pattern=[&](std::size_t r){for(std::size_t i=0;i<n;++i)if(f[start+i]!=(i<r))return false;return true;};
    if(n>=8&&g4&&pattern(4)){return "gpc4";}
    if(n>=4&&g2&&pattern(2)){return "gpc2";}
    if(n>1&&spc&&pattern(1)){return "spc";}
    if(n==1){return f[start]?"frozen":"information";}
    return {};
}
void leaf_range(const std::vector<double>&llrs,std::vector<DecoderPath>&paths,const std::vector<bool>&flags,std::size_t start,std::size_t size,std::size_t L,DecoderStats&stats){
    for(std::size_t i=start;i<start+size;++i){std::vector<DecoderPath> expanded;++stats.leaf_steps;
        for(const auto&path:paths){auto a=polar_scl::node_alpha(llrs,path.u,i,1)[0];for(Bit b=0;b<(flags[i]?1:2);++b){auto p=path;p.metric+=(b==static_cast<Bit>(a<0)?0:std::abs(a));p.u.push_back(b);expanded.push_back(p);}if(!flags[i])++stats.path_splits;}
        paths=polar_scl::prune_paths(std::move(expanded),L);if(!flags[i])++stats.sorts;}
}
}
void NodeCostEstimate::validate()const{if(size==0||list_size==0||!std::isfinite(mean_work)||!std::isfinite(p99_work)||mean_work<0||p99_work<mean_work||!std::isfinite(fallback_rate)||fallback_rate<0||fallback_rate>1)throw std::invalid_argument("invalid node cost estimate");}
void ScheduleCostModel::validate()const{if(!std::isfinite(tail_weight)||!std::isfinite(fallback_penalty)||tail_weight<0||fallback_penalty<0)throw std::invalid_argument("invalid risk weights");std::set<std::tuple<std::string,std::size_t,std::size_t>>keys;for(const auto&p:profiles){p.validate();if(!keys.emplace(p.kind,p.size,p.list_size).second)throw std::invalid_argument("duplicate cost profile");}}
std::optional<double> ScheduleCostModel::score(const std::string&kind,std::size_t n,std::size_t L)const{for(const auto&p:profiles)if(p.kind==kind&&p.size==n&&p.list_size==L)return p.mean_work+tail_weight*(p.p99_work-p.mean_work)+fallback_penalty*p.fallback_rate;return {};}
std::vector<DecodeOp> rate1_schedule(const std::vector<bool>&f,bool spc,bool g2,bool g4){validate_flags(f);std::vector<DecodeOp>out;
    std::function<void(std::size_t,std::size_t)>visit=[&](auto start,auto n){auto kind=kind_at(f,start,n,spc,g2,g4);if(!kind.empty()){out.push_back({kind,start,n});return;}visit(start,n/2);visit(start+n/2,n/2);};visit(0,f.size());return out;}
ScheduleDecision risk_aware_schedule(const std::vector<bool>&f,std::size_t L,const ScheduleCostModel&cost,bool spc,bool g2,bool g4){validate_flags(f);cost.validate();if(L==0)throw std::invalid_argument("list_size must be positive");
    std::function<ScheduleDecision(std::size_t,std::size_t)>solve=[&](auto start,auto n)->ScheduleDecision{if(n==1)return {{{f[start]?"frozen":"information",start,1}},f[start]?cost.frozen_leaf_cost:cost.information_leaf_cost,{}};
        auto left=solve(start,n/2),right=solve(start+n/2,n/2);left.operations.insert(left.operations.end(),right.operations.begin(),right.operations.end());left.estimated_cost+=right.estimated_cost;left.rejected_constituents.insert(left.rejected_constituents.end(),right.rejected_constituents.begin(),right.rejected_constituents.end());
        auto kind=kind_at(f,start,n,spc,g2,g4);if(kind.empty())return left;auto score=cost.score(kind,n,L);if(score&&*score<=left.estimated_cost)return {{{kind,start,n}},*score,{}};
        left.rejected_constituents.push_back(kind+"@"+std::to_string(start)+"+"+std::to_string(n)+(score?":cost_exceeds_children":":unprofiled"));return left;};return solve(0,f.size());}
void DecoderStats::include(const std::string&kind,const SearchStats&s){
    if(kind=="rate1"){++rate1_nodes;rate1_candidates+=s.candidates_materialized;rate1_frontier_expansions+=s.frontier_expansions;rate1_fallbacks+=s.fallback_used;}
    else if(kind=="spc"){++spc_nodes;spc_candidates+=s.candidates_materialized;spc_frontier_expansions+=s.frontier_expansions;spc_fallbacks+=s.fallback_used;spc_parity_rejections+=s.parity_rejections;}
    else if(kind=="gpc2"){++gpc2_nodes;gpc2_candidates+=s.candidates_materialized;gpc2_frontier_expansions+=s.frontier_expansions;gpc2_fallbacks+=s.fallback_used;gpc2_syndrome_rejections+=s.syndrome_rejections;}
    else if(kind=="gpc4"){++gpc4_nodes;gpc4_candidates+=s.candidates_materialized;gpc4_frontier_expansions+=s.frontier_expansions;gpc4_fallbacks+=s.fallback_used;gpc4_syndrome_rejections+=s.syndrome_rejections;}}
DecodeResult rate1_scl_decode(const std::vector<double>&llrs,const std::vector<bool>&flags,std::size_t L,const Rate1Policy&policy,bool spc,bool g2,bool g4,const PolarDecoderIR*supplied,const ScheduleCostModel*cost){
    polar_scl::validate_decoder_inputs(llrs,flags,L);auto ir=supplied?*supplied:default_polar_decoder_ir(policy);validate_polar_decoder_ir(ir);spc=spc&&ir.enable_spc;g2=g2&&ir.enable_gpc2;g4=g4&&ir.enable_gpc4;
    auto operations=cost?risk_aware_schedule(flags,L,*cost,spc,g2,g4).operations:rate1_schedule(flags,spc,g2,g4);DecodeResult result;result.paths.push_back({0,{}});auto&paths=result.paths;auto&stats=result.stats;
    for(const auto&op:operations){for(const auto&path:paths)if(path.u.size()!=op.start)throw std::logic_error("nonsequential decode schedule");
        if(op.kind=="rate0"){++stats.rate0_nodes;for(auto&path:paths){auto alpha=polar_scl::node_alpha(llrs,path.u,op.start,op.size);std::vector<double>terms;for(auto value:alpha)if(value<0)terms.push_back(std::abs(value));path.metric+=metric_sum(terms);path.u.insert(path.u.end(),op.size,0);}paths=polar_scl::prune_paths(std::move(paths),L);continue;}
        if(op.kind=="frozen"||op.kind=="information"){leaf_range(llrs,paths,flags,op.start,1,L,stats);continue;}
        std::vector<ParentPath>parents;for(std::size_t i=0;i<paths.size();++i)parents.push_back({i,paths[i].metric,polar_scl::node_alpha(llrs,paths[i].u,op.start,op.size)});
        auto rank=[&](const Rate1Candidate&c){auto u=polar_transform(c.bits);std::vector<std::uint64_t>key(paths[c.parent_id].u.begin(),paths[c.parent_id].u.end());key.insert(key.end(),u.begin(),u.end());return polar_scl::GeneralRankKey{c.path_metric,key};};
        auto hidden=[&](const ParentPath&p,double metric){std::vector<std::uint64_t>key(paths[p.parent_id].u.begin(),paths[p.parent_id].u.end());key.insert(key.end(),op.size,0);return polar_scl::GeneralRankKey{metric,key};};
        const auto&algorithm=op.kind=="rate1"?ir.rate1:op.kind=="spc"?ir.spc:op.kind=="gpc2"?ir.gpc2:ir.gpc4;
        auto materialized=execute_algorithm_ir_integrated(algorithm,parents,L,rank,hidden);stats.include(op.kind,materialized.stats);
        if(!materialized.candidates||materialized.stats.fallback_used){leaf_range(llrs,paths,flags,op.start,op.size,L,stats);continue;}
        std::vector<DecoderPath>expanded;for(const auto&c:*materialized.candidates){auto u=polar_transform(c.bits);auto checks=op.kind=="spc"?1u:op.kind=="gpc2"?2u:op.kind=="gpc4"?4u:0u;
            if(std::any_of(u.begin(),u.begin()+checks,[](Bit b){return b!=0;})){throw std::logic_error("constituent violated frozen bits");}
            auto path=paths[c.parent_id];path.metric=c.path_metric;path.u.insert(path.u.end(),u.begin(),u.end());expanded.push_back(std::move(path));}
        paths=polar_scl::prune_paths(std::move(expanded),L);++stats.sorts;}
    return result;
}
}
