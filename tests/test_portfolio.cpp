#include "test_util.hpp"
#include <random>
using namespace polar_scl;
int main()
{
    return run_test_main("portfolio",[]
    {
        PortfolioEvaluationConfig ec;
        ec.node_sizes={4,8};
        ec.list_sizes={2,4};
        ec.trials_per_case=8;
        ec.quantization_step=.25;
        ec.seed=3;
        auto r=evaluate_selector(SelectorGenome{},ec);
        CHECK_TRUE(r.is_exact());
        CHECK_TRUE(!r.engine_counts.empty());
        std::mt19937_64 rng(53);
        std::normal_distribution<double>nd(0,3);
        PolarAlgorithmIR alg;
        for(int t=0;t<15;++t)
        {
            auto f=pw_frozen_pattern(32,20);
            std::vector<double>llrs(32);
            for(double&v:llrs)v=nd(rng);
            auto ref=scl_decode(llrs,f,4);
            auto a=portfolio_scl_decode(llrs,f,4,alg);
            check_paths_equal(a.paths,ref.paths,1e-11);
        }

        SelectorGenome g;
        g.small_node_threshold=8;
        g.small_engine=Rate1Engine::FastSSCLBounded;
        g.easy_engine=Rate1Engine::LazyExact;
        g.hard_engine=Rate1Engine::FastSSCLBounded;
        alg.selector=g;
        std::uniform_int_distribution<int>q(-5,5);
        for(int t=0;t<15;++t)
        {
            auto f=pw_frozen_pattern(16,10);
            std::vector<double>llrs(16);
            for(double&v:llrs)v=q(rng);
            auto ref=scl_decode(llrs,f,4);
            auto a=portfolio_scl_decode(llrs,f,4,alg);
            check_paths_equal(a.paths,ref.paths,0.0);
        }

        PortfolioEvaluationConfig e;
        e.node_sizes={4,8};
        e.list_sizes={2,4};
        e.trials_per_case=3;
        e.quantization_step=.25;
        e.seed=11;
        SelectorEvolutionConfig evo;
        evo.population_size=4;
        evo.offspring_per_generation=3;
        evo.generations=2;
        evo.seed=13;
        auto a=evolutionary_selector_search({},evo,e);
        auto bres=evolutionary_selector_search({},evo,e);
        CHECK_EQ(a.evaluated.size(),bres.evaluated.size());
        CHECK_TRUE(a.best.is_exact());
        CHECK_TRUE(close_metric(a.best.mean_proxy_cost,bres.best.mean_proxy_cost,0.0));
        CHECK_TRUE(a.best.algorithm.selector==bres.best.algorithm.selector);
        for(const auto&x:a.evaluated)CHECK_TRUE(x.is_exact());
    }

    );
}
