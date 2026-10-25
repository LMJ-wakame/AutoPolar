#include "polar_scl/polar_scl.hpp"
#include <iomanip>
#include <iostream>
using namespace polar_scl;
int main()
{
    SelectorEvolutionConfig evo;
    evo.population_size=10;
    evo.offspring_per_generation=10;
    evo.generations=6;
    evo.seed=17;
    PortfolioEvaluationConfig eval;
    eval.node_sizes={4,8,16,32};
    eval.list_sizes={2,4,8};
    eval.trials_per_case=15;
    eval.quantization_step=.25;
    eval.seed=42;
    auto r=evolutionary_selector_search({},evo,eval);
    std::cout<<"Evaluated "<<r.evaluated.size()<<" selector genomes\nBest proxy cost by generation:";
    for(double x:r.best_cost_by_generation)std::cout<<' '<<std::fixed<<std::setprecision(3)<<x;
    std::cout<<"\nBest exact selector: mean_proxy="<<r.best.mean_proxy_cost<<" p99="<<r.best.p99_proxy_cost<<" mean_candidates="<<r.best.mean_candidates<<" fallback="<<r.best.fallback_rate<<"\n";
    const auto&g=r.best.algorithm.selector;
    std::cout<<"small<= "<<g.small_node_threshold<<" -> "<<to_string(g.small_engine)<<", easy -> "<<to_string(g.easy_engine)<<", hard -> "<<to_string(g.hard_engine)<<", gap>= "<<g.reliability_gap_threshold<<", PMgap>= "<<g.parent_metric_gap_threshold<<"\n";
}
