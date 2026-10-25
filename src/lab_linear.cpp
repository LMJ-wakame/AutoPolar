#include "polar_scl_lab/linear.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <stdexcept>

namespace polar_scl_lab {
namespace {
using Mask=std::uint64_t;
constexpr double inf=std::numeric_limits<double>::infinity();
using Row=std::vector<double>;
// Bounds round down, while emitted candidate metrics use accurate sums.
double lower_add(double a,double b){if(!std::isfinite(a)||!std::isfinite(b))return a+b;return std::nextafter(a+b,-inf);}
double lower_nonnegative_add(double a,double b){return std::max(0.0,lower_add(a,b));}
double accurate_sum(const Row& values) {
    Row partials;
    for(double x: values) {
        std::size_t i=0;
        for(double y: partials) {
            if(std::abs(x)<std::abs(y)) std::swap(x,y);
            const double hi=x+y, lo=y-(hi-x);
            if(lo!=0.0) partials[i++]=lo;
            x=hi;
        }
        partials.resize(i); partials.push_back(x);
    }
    double hi=0.0,lo=0.0;
    while(!partials.empty()) {
        const double x=hi,y=partials.back(); partials.pop_back();
        hi=x+y; lo=y-(hi-x); if(lo!=0.0) break;
    }
    if(!partials.empty() && ((lo<0 && partials.back()<0)||(lo>0 && partials.back()>0))) {
        const double y=lo*2.0,x=hi+y; if(x-hi==y) hi=x;
    }
    return hi;
}
bool less_candidate(const Rate1Candidate&a,const Rate1Candidate&b) { return polar_scl::local_rank_key(a)<polar_scl::local_rank_key(b); }
void trim(std::vector<Rate1Candidate>& v,std::size_t L) { std::sort(v.begin(),v.end(),less_candidate); if(v.size()>L)v.resize(L); }
void validate(const std::vector<ParentPath>& p,std::size_t L,const LinearConstraintIR& h) {
    h.validate(); polar_scl::validate_parent_paths(p,L);
    if(p.front().llrs.size()!=h.num_variables) throw std::invalid_argument("constraint width differs from parent size");
}
Mask syndrome(Mask m,const std::vector<Mask>& cols) { Mask s=0;for(std::size_t i=0;i<cols.size();++i)if((m>>i)&1ULL)s^=cols[i];return s; }
Mask target(const ParentPath&p,const std::vector<Mask>&cols) { Mask s=0;for(std::size_t i=0;i<cols.size();++i)if(p.llrs[i]<0)s^=cols[i];return s; }
std::size_t power(std::size_t k) { if(k>=63)throw std::length_error("enumerated frontier exceeds addressable mask range");return std::size_t{1}<<k; }
void merge(SearchStats&a,const SearchStats&b) {
    a.candidates_materialized+=b.candidates_materialized;a.subset_states_pushed+=b.subset_states_pushed;
    a.subset_states_popped+=b.subset_states_popped;a.global_heap_pushes+=b.global_heap_pushes;a.global_heap_pops+=b.global_heap_pops;
    a.syndrome_rejections+=b.syndrome_rejections;a.metric_evaluations+=b.metric_evaluations;a.constraint_checks+=b.constraint_checks;
    a.feasible_candidates+=b.feasible_candidates;a.dp_state_relaxations+=b.dp_state_relaxations;
}
std::vector<std::size_t> positions(const ParentPath&p) {
    std::vector<std::size_t> pos(p.llrs.size());std::iota(pos.begin(),pos.end(),0);
    std::sort(pos.begin(),pos.end(),[&](auto a,auto b){return std::pair{std::abs(p.llrs[a]),a}<std::pair{std::abs(p.llrs[b]),b};});return pos;
}
Mask original(Mask sorted,const std::vector<std::size_t>&pos) { Mask out=0;for(std::size_t i=0;i<pos.size();++i)if((sorted>>i)&1ULL)out|=1ULL<<pos[i];return out; }
struct SubsetStream {
    const ParentPath* parent;SearchStats* stats;std::vector<std::size_t> pos;Row weights;
    using Entry=std::tuple<double,Mask,std::size_t,Mask>;
    std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> heap;
    bool empty_yielded=false,initialized=false;
    SubsetStream(const ParentPath&p,SearchStats&s):parent(&p),stats(&s),pos(positions(p)){for(auto i:pos)weights.push_back(std::abs(p.llrs[i]));}
    void push(double penalty,Mask m,std::size_t last){heap.emplace(penalty,original(m,pos),last,m);++stats->subset_states_pushed;}
    std::optional<Rate1Candidate> next(){
        if(!empty_yielded){empty_yielded=true;return make_candidate(*parent,0);}
        if(!initialized){initialized=true;push(weights[0],1,0);}
        if(heap.empty())return {};
        auto [penalty,m,last,sorted]=heap.top();heap.pop();++stats->subset_states_popped;
        auto next=last+1;if(next<pos.size()){auto bit=1ULL<<next;push(penalty+weights[next],sorted|bit,next);push(penalty-weights[last]+weights[next],(sorted^(1ULL<<last))|bit,next);}
        return make_candidate(*parent,m);
    }
};
template<class Stream> Rate1SearchResult merge_streams(std::vector<Stream>&streams,std::size_t L,SearchStats&stats) {
    struct Entry {Rate1Candidate candidate;std::size_t stream;};
    auto greater=[](const Entry&a,const Entry&b){return less_candidate(b.candidate,a.candidate);};
    std::priority_queue<Entry,std::vector<Entry>,decltype(greater)> heap(greater);
    for(std::size_t i=0;i<streams.size();++i)if(auto c=streams[i].next()){heap.push({*c,i});++stats.global_heap_pushes;}
    std::vector<Rate1Candidate> result;
    while(!heap.empty()&&result.size()<L){auto e=heap.top();heap.pop();++stats.global_heap_pops;result.push_back(e.candidate);
        if(auto c=streams[e.stream].next()){heap.push({*c,e.stream});++stats.global_heap_pushes;}}
    trim(result,L);return {result,stats};
}
Row advance(const Row&prev,double weight,Mask column){Row next(prev.size());for(std::size_t s=0;s<prev.size();++s)next[s]=std::min(prev[s],lower_nonnegative_add(prev[s^column],weight));return next;}
struct Frontier {
    const ParentPath* parent;std::vector<std::size_t> pos;Row weights;std::vector<Mask> cols;Mask required;
    std::vector<Row> prefix,hidden;std::vector<std::vector<std::size_t>> counts;Row penalties;std::size_t active;
    Frontier(const ParentPath&p,const LinearConstraintIR&h,std::size_t initial,bool parity_only=false):parent(&p),pos(positions(p)),active(std::min(initial,pos.size())){
        auto raw=h.columns();required=target(p,raw);for(auto i:pos){weights.push_back(std::abs(p.llrs[i]));cols.push_back(raw[i]);}
        if(parity_only){
            counts.resize(pos.size()+1,std::vector<std::size_t>(2));counts[0][0]=1;
            for(std::size_t k=1;k<=pos.size();++k)counts[k][0]=counts[k][1]=std::size_t{1}<<(k-1);
            for(std::size_t k=0;k<=pos.size();++k){double best=inf;
                if(k<pos.size()){
                    if(required==1)best=weights[k];
                    else if(k>0)best=lower_nonnegative_add(weights[k],weights[0]);
                    if(k+1<pos.size()){
                        if(required==0)best=std::min(best,lower_nonnegative_add(weights[k],weights[k+1]));
                        else if(k>0)best=std::min(best,lower_nonnegative_add(lower_nonnegative_add(weights[k],weights[k+1]),weights[0]));
                    }
                }
                penalties.push_back(best);
            }
            return;
        }
        if(h.rank()>16)throw std::length_error("frontier DP rank limit is 16");
        const auto S=h.state_count(),n=pos.size();Row zero(S,inf);zero[0]=0;prefix.push_back(zero);counts.emplace_back(S,0);counts[0][0]=1;
        for(std::size_t i=0;i<n;++i){prefix.push_back(advance(prefix.back(),weights[i],cols[i]));std::vector<std::size_t> row(S);
            for(std::size_t s=0;s<S;++s){auto a=counts.back()[s],b=counts.back()[s^cols[i]];row[s]=a>SIZE_MAX-b?SIZE_MAX:a+b;}counts.push_back(row);}
        std::vector<Row> any(n+1);hidden.resize(n+1);any[n]=zero;hidden[n]=Row(S,inf);
        for(std::size_t i=n;i-->0;){any[i]=advance(any[i+1],weights[i],cols[i]);hidden[i].resize(S);
            for(std::size_t s=0;s<S;++s)hidden[i][s]=std::min(hidden[i+1][s],lower_nonnegative_add(any[i+1][s^cols[i]],weights[i]));}
        for(std::size_t k=0;k<=n;++k){double best=inf;for(std::size_t s=0;s<S;++s)best=std::min(best,lower_nonnegative_add(hidden[k][s],prefix[k][required^s]));penalties.push_back(best);}
    }
    std::size_t count(std::size_t k)const{return counts[k][required];}
    std::optional<double> lower()const{if(active==pos.size()||!std::isfinite(penalties[active]))return {};return lower_add(parent->path_metric,penalties[active]);}
};
struct BudgetExceeded {};
struct AStar {
    const ParentPath* parent;SearchStats* stats;std::optional<std::size_t> budget;
    std::vector<std::size_t> pos;Row weights;std::vector<Mask> cols;Mask required;std::vector<Row> suffix;
    using Entry=std::tuple<double,double,std::size_t,Mask,Mask>;
    std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> heap;
    struct Greater {bool operator()(const Rate1Candidate&a,const Rate1Candidate&b)const{return less_candidate(b,a);}};
    std::priority_queue<Rate1Candidate,std::vector<Rate1Candidate>,Greater> pending;
    AStar(const ParentPath&p,const LinearConstraintIR&h,SearchStats&s,std::optional<std::size_t>b):parent(&p),stats(&s),budget(b),pos(positions(p)){
        auto raw=h.columns();required=target(p,raw);for(auto i:pos){weights.push_back(std::abs(p.llrs[i]));cols.push_back(raw[i]);}
        if(h.rank()>16)throw std::length_error("A* DP rank limit is 16");
        auto S=h.state_count();suffix.resize(pos.size()+1);suffix.back()=Row(S,inf);suffix.back()[0]=0;
        for(std::size_t i=pos.size();i-->0;)suffix[i]=advance(suffix[i+1],weights[i],cols[i]);
        stats->dp_state_relaxations+=pos.size()*S;
        if(std::isfinite(suffix[0][required])){heap.emplace(suffix[0][required],0,0,0,0);++stats->subset_states_pushed;}
    }
    void push(double penalty,std::size_t d,Mask s,Mask m){double heuristic=suffix[d][required^s];if(!std::isfinite(heuristic))return;
        if(budget&&stats->subset_states_pushed+1>*budget)throw BudgetExceeded{};
        heap.emplace(lower_nonnegative_add(penalty,heuristic),penalty,d,s,m);++stats->subset_states_pushed;}
    std::optional<Rate1Candidate> next(){while(true){double cutoff=pending.empty()?inf:pending.top().path_metric;
        if(!pending.empty()&&(heap.empty()||lower_add(parent->path_metric,std::get<0>(heap.top()))>cutoff)){auto c=pending.top();pending.pop();return c;}
        if(heap.empty())return {};
        auto [bound,penalty,d,s,m]=heap.top();(void)bound;heap.pop();++stats->subset_states_popped;
        if(d==pos.size()){if(s!=required)throw std::logic_error("infeasible A* terminal");auto c=make_candidate(*parent,original(m,pos));
            ++stats->candidates_materialized;++stats->metric_evaluations;++stats->constraint_checks;++stats->feasible_candidates;pending.push(c);continue;}
        push(penalty,d+1,s,m);push(lower_nonnegative_add(penalty,weights[d]),d+1,s^cols[d],m|(1ULL<<d));}}
};
}
Rate1Candidate make_candidate(const ParentPath&p,std::uint64_t m){Rate1Candidate c;c.parent_id=p.parent_id;c.flip_mask=m;Row terms;
    for(std::size_t i=0;i<p.llrs.size();++i){auto flip=static_cast<Bit>((m>>i)&1ULL);c.bits.push_back(static_cast<Bit>((p.llrs[i]<0)^flip));if(flip)terms.push_back(std::abs(p.llrs[i]));}
    c.penalty=accurate_sum(terms);c.path_metric=p.path_metric+c.penalty;return c;}
double metric_sum(const std::vector<double>&values){return accurate_sum(values);}
LinearConstraintIR gpc_constraint(std::size_t n,std::size_t r){if(r<2||(r&(r-1))||n<=r||(n&(n-1))||n>64)throw std::invalid_argument("invalid GPC dimensions");
    LinearConstraintIR h{"gpc"+std::to_string(r),n,{}};for(std::size_t j=0;j<r;++j){Mask row=0;for(std::size_t i=j;i<n;i+=r)row|=1ULL<<i;h.parity_checks.push_back(row);}h.validate();return h;}
ConstraintProofMetadata proof_metadata(const LinearConstraintIR&h){return {h.name,"D_next[s]=min(D[s],D[s xor column]+weight)",h.rank(),h.state_count(),
    {"compiled columns preserve parity equations","visible candidates satisfy target syndrome","hidden excludes empty subset","bound covers every unseen feasible candidate","strict total-order termination","budget invokes exact fallback"}};}
Rate1SearchResult exhaustive_linear_top_l(const std::vector<ParentPath>&p,std::size_t L,const LinearConstraintIR&h){validate(p,L,h);if(h.num_variables>24)throw std::length_error("exhaustive oracle limited to 24 bits");Rate1SearchResult out;out.stats.algorithm="linear_exhaustive:"+h.name;auto cols=h.columns();
    for(const auto&parent:p){auto required=target(parent,cols);for(Mask m=0;m<power(h.num_variables);++m){++out.stats.constraint_checks;if(syndrome(m,cols)==required){out.candidates.push_back(make_candidate(parent,m));++out.stats.metric_evaluations;++out.stats.feasible_candidates;}}}
    out.stats.candidates_materialized=out.candidates.size();out.stats.syndrome_rejections=out.stats.constraint_checks-out.candidates.size();trim(out.candidates,L);return out;}
Rate1SearchResult lazy_exact_linear_top_l(const std::vector<ParentPath>&p,std::size_t L,const LinearConstraintIR&h){validate(p,L,h);SearchStats stats;stats.algorithm="linear_lazy_exact:"+h.name;
    struct Stream {SubsetStream subset;std::vector<Mask> cols;Mask required;SearchStats* stats;
        Stream(const ParentPath&p,const LinearConstraintIR&h,SearchStats&s):subset(p,s),cols(h.columns()),required(target(p,cols)),stats(&s){}
        std::optional<Rate1Candidate> next(){while(auto c=subset.next()){++stats->candidates_materialized;++stats->metric_evaluations;++stats->constraint_checks;
            if(syndrome(c->flip_mask,cols)==required){++stats->feasible_candidates;return c;}++stats->syndrome_rejections;}return {};}};
    std::vector<Stream> streams;for(const auto&parent:p)streams.emplace_back(parent,h,stats);return merge_streams(streams,L,stats);}
Rate1SearchResult astar_linear_top_l(const std::vector<ParentPath>&p,std::size_t L,const LinearConstraintIR&h,std::optional<std::size_t> budget){validate(p,L,h);if(budget&&*budget==0)throw std::invalid_argument("state_budget must be positive");SearchStats stats;stats.algorithm="linear_astar_exact:"+h.name;
    try{std::vector<AStar> streams;for(const auto&parent:p)streams.emplace_back(parent,h,stats,budget);return merge_streams(streams,L,stats);}
    catch(const BudgetExceeded&){auto exact=lazy_exact_linear_top_l(p,L,h);stats.fallback_used=true;stats.fallback_reason="A* state budget exceeded";merge(stats,exact.stats);exact.stats=stats;return exact;}}
static Rate1SearchResult adaptive_impl(const std::vector<ParentPath>&p,std::size_t L,const LinearConstraintIR&h,const Rate1Policy&policy,CandidateRankFn rank,HiddenRankFn hidden_rank,bool parity_only){
    validate(p,L,h);validate_policy(policy);if(!rank)rank=polar_scl::general_local_rank_key;
    SearchStats stats;stats.algorithm="linear_adaptive_certified:"+h.name;std::vector<Frontier> states;
    for(const auto&parent:p){states.emplace_back(parent,h,policy.initial_active_bits,parity_only);if(!parity_only)stats.dp_state_relaxations+=(5*h.num_variables+1)*h.state_count();}
    auto lower=[&](const Frontier&s)->std::optional<polar_scl::GeneralRankKey>{auto value=s.lower();if(!value)return {};if(hidden_rank)return hidden_rank(*s.parent,*value);return polar_scl::GeneralRankKey{*value,{s.parent->parent_id,0}};};
    auto fallback=[&](const char*reason){auto exact=lazy_exact_linear_top_l(p,L,h);stats.fallback_used=true;stats.fallback_reason=reason;merge(stats,exact.stats);exact.stats=stats;return exact;};
    std::vector<Rate1Candidate> visible;std::size_t initial=0;for(const auto&s:states){auto c=s.count(s.active);if(c>policy.candidate_budget-initial)return fallback("initial linear frontier exceeds candidate budget");initial+=c;}
    auto emit=[&](Frontier&s,std::size_t begin){auto end=power(s.active);for(std::size_t m=begin;m<end;++m){++stats.constraint_checks;
        if(syndrome(m,s.cols)!=s.required){continue;}
        visible.push_back(make_candidate(*s.parent,original(m,s.pos)));++stats.metric_evaluations;++stats.feasible_candidates;++stats.candidates_materialized;}
        stats.per_parent_active_bits[s.parent->parent_id]=s.active;};
    for(auto&s:states)emit(s,0);
    // Nullspace count is needed when the requested list exceeds all feasible words.
    auto freedom=h.num_variables-h.rank();std::size_t target_count=L;
    if(freedom<63&&p.size()<=SIZE_MAX/power(freedom))target_count=std::min(L,p.size()*power(freedom));
    while(true){std::sort(visible.begin(),visible.end(),[&](const auto&a,const auto&b){return rank(a)<rank(b);});++stats.certificate_checks;
        std::optional<polar_scl::GeneralRankKey> best;for(const auto&s:states)if(auto key=lower(s);key&&(!best||*key<*best))best=key;
        const bool enough=visible.size()>=target_count;
        if(enough){auto worst=rank(visible[target_count-1]);stats.certificate_worst_key=worst;stats.certificate_hidden_key=best;
            if(!best||*best>worst){visible.resize(target_count);return {visible,stats};}}
        std::vector<std::size_t> competitive;for(std::size_t i=0;i<states.size();++i)if(auto key=lower(states[i]);key&&(!enough||*key<=rank(visible[target_count-1])))competitive.push_back(i);
        if(competitive.empty())return fallback("no expandable linear frontier before certificate");
        if(policy.expansion_mode==ExpansionMode::MinimumBound){auto chosen=*std::min_element(competitive.begin(),competitive.end(),[&](auto a,auto b){return *lower(states[a])<*lower(states[b]);});competitive={chosen};}
        std::size_t additional=0;for(auto i:competitive){auto&s=states[i];auto next=std::min(s.active+policy.growth_step,s.pos.size());auto count=s.count(next)-s.count(s.active);
            if(count>policy.candidate_budget-std::min(policy.candidate_budget,stats.candidates_materialized+additional)){return fallback("linear frontier expansion exceeds candidate budget");}
            additional+=count;}
        for(auto i:competitive){auto&s=states[i];auto old=s.active;s.active=std::min(s.active+policy.growth_step,s.pos.size());emit(s,power(old));++stats.frontier_expansions;}
    }
}
Rate1SearchResult adaptive_linear_top_l(const std::vector<ParentPath>&p,std::size_t L,const LinearConstraintIR&h,const Rate1Policy&policy,CandidateRankFn rank,HiddenRankFn hidden){return adaptive_impl(p,L,h,policy,rank,hidden,false);}
static LinearConstraintIR spc_h(const std::vector<ParentPath>&p){if(p.empty())throw std::invalid_argument("empty parents");auto n=p.front().llrs.size();if(n==0||n>64)throw std::invalid_argument("SPC width");return {"spc",n,{n==64?~0ULL:(1ULL<<n)-1}};}
Rate1SearchResult exhaustive_spc_top_l(const std::vector<ParentPath>&p,std::size_t L){auto out=exhaustive_linear_top_l(p,L,spc_h(p));out.stats.parity_rejections=out.stats.syndrome_rejections;return out;}
Rate1SearchResult lazy_exact_spc_top_l(const std::vector<ParentPath>&p,std::size_t L){auto out=lazy_exact_linear_top_l(p,L,spc_h(p));out.stats.parity_rejections=out.stats.syndrome_rejections;return out;}
Rate1SearchResult adaptive_spc_top_l(const std::vector<ParentPath>&p,std::size_t L,const Rate1Policy&policy,CandidateRankFn rank,HiddenRankFn hidden){auto out=adaptive_impl(p,L,spc_h(p),policy,rank,hidden,true);out.stats.parity_rejections=out.stats.syndrome_rejections;return out;}
Rate1SearchResult exhaustive_gpc_top_l(const std::vector<ParentPath>&p,std::size_t L,std::size_t r){if(p.empty())throw std::invalid_argument("empty parents");return exhaustive_linear_top_l(p,L,gpc_constraint(p.front().llrs.size(),r));}
Rate1SearchResult lazy_exact_gpc_top_l(const std::vector<ParentPath>&p,std::size_t L,std::size_t r){if(p.empty())throw std::invalid_argument("empty parents");return lazy_exact_linear_top_l(p,L,gpc_constraint(p.front().llrs.size(),r));}
Rate1SearchResult adaptive_gpc_top_l(const std::vector<ParentPath>&p,std::size_t L,std::size_t r,const Rate1Policy&policy,CandidateRankFn rank,HiddenRankFn hidden){if(p.empty())throw std::invalid_argument("empty parents");return adaptive_linear_top_l(p,L,gpc_constraint(p.front().llrs.size(),r),policy,rank,hidden);}
Rate1SearchResult CompiledLinearAlgorithm::execute(const std::vector<ParentPath>&p,std::size_t L)const{switch(engine){case LinearEngine::LazyExact:return lazy_exact_linear_top_l(p,L,constraint);case LinearEngine::CertifiedFrontier:return adaptive_linear_top_l(p,L,constraint,policy);case LinearEngine::SyndromeAStar:return astar_linear_top_l(p,L,constraint,state_budget);case LinearEngine::SparseTrellis:{auto ir=search_trellis_ir(constraint,128,73);return linear_trellis_top_l(p,L,ir);}}throw std::logic_error("unknown linear engine");}
CompiledLinearAlgorithm::CompiledLinearAlgorithm(LinearConstraintIR h,LinearEngine e,std::optional<Rate1Policy>p,std::optional<std::size_t>b):constraint(std::move(h)),engine(e),policy(p.value_or(Rate1Policy{})),state_budget(b){constraint.validate();validate_policy(policy);if(p&&engine!=LinearEngine::CertifiedFrontier)throw std::invalid_argument("policy applies only to certified frontier");if(b&&(*b==0||engine!=LinearEngine::SyndromeAStar))throw std::invalid_argument("state budget applies to A* and must be positive");}
CompiledLinearAlgorithm compile_linear_algorithm_ir(const LinearConstraintIR&h,LinearEngine engine,std::optional<Rate1Policy>policy,std::optional<std::size_t>budget){return CompiledLinearAlgorithm(h,engine,policy,budget);}
}
