#include "polar_scl/polar_scl.hpp"
#include <iomanip>
#include <iostream>
using namespace polar_scl;
static void print(const char*name,const PortfolioResult&r)
{
    std::cout<<std::left<<std::setw(22)<<name<<std::right<<" exact="<<r.exact_cases<<'/'<<r.total_cases<<" mean_proxy="<<std::fixed<<std::setprecision(3)<<r.mean_proxy_cost<<" p99="<<r.p99_proxy_cost<<" mean_candidates="<<r.mean_candidates<<" fallback="<<r.fallback_rate<<" engines=";
    for(auto[e,c]:r.engine_counts)std::cout<<to_string(e)<<':'<<c<<' ';
    std::cout<<'\n';
}

int main()
{
    PortfolioEvaluationConfig cfg;
    cfg.node_sizes={4,8,16,32};
    cfg.list_sizes={2,4,8};
    cfg.trials_per_case=25;
    cfg.quantization_step=.25;
    cfg.seed=42;
    for(auto e:{Rate1Engine::LazyExact,Rate1Engine::FastSSCLBounded,Rate1Engine::CertifiedFrontier})
    {
        PolarAlgorithmIR a;
        a.selector_mode=SelectorMode::Fixed;
        a.fixed_rate1.engine=e;
        if(e==Rate1Engine::CertifiedFrontier)
        {
            a.fixed_rate1.initial_active_bits=0;
            a.fixed_rate1.candidate_budget=256;
        }

        auto r=evaluate_algorithm(a,cfg);
        print(to_string(e).c_str(),r);
    }

    PolarAlgorithmIR dynamic;
    auto r=evaluate_algorithm(dynamic,cfg);
    print("default_rule_selector",r);
}
