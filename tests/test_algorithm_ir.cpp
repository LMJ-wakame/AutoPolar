#include "test_util.hpp"
#include <random>
#include <set>
using namespace polar_scl;
int main()
{
    return run_test_main("algorithm_ir",[]
    {
        std::mt19937_64 rng(2026);
        std::uniform_int_distribution<int>v(-4,4),m(0,3);
        std::vector<Rate1NodeIR> engines{{Rate1Engine::Exhaustive},{Rate1Engine::LazyExact},{Rate1Engine::FastSSCLBounded},{Rate1Engine::CertifiedFrontier,0,1,ExpansionMode::MinimumBound,256}};
        for(std::size_t ns:{1u,2u,4u,6u})for(std::size_t L:
        {
            1u,2u,4u
        }

        )for(int t=0;t<20;++t)
        {
            std::vector<ParentPath>p;
            for(std::size_t id=0;id<L;++id)
            {
                std::vector<double>x(ns);
                for(double&a:x)a=v(rng);
                p.push_back(
                {
                    id,static_cast<double>(m(rng)),x
                }

                );
            }

            auto ex=exhaustive_rate1_top_l(p,L);
            for(auto ir:engines)
            {
                auto a=execute_rate1_ir(ir,p,L);
                check_candidates_equal(a.candidates,ex.candidates);
            }
        }

        SelectorGenome g;
        g.small_node_threshold=2;
        g.reliability_gap_threshold=.5;
        g.parent_metric_gap_threshold=.25;
        g.small_engine=Rate1Engine::Exhaustive;
        g.easy_engine=Rate1Engine::CertifiedFrontier;
        g.hard_engine=Rate1Engine::LazyExact;
        PolarAlgorithmIR alg;
        alg.selector=g;
        std::set<Rate1Engine>used;
        std::uniform_real_distribution<double>ur(-4,4),pm(0,1);
        for(int t=0;t<60;++t)
        {
            std::vector<ParentPath>p;
            for(std::size_t id=0;id<4;++id)
            {
                std::vector<double>x(8);
                for(double&a:x)a=ur(rng);
                p.push_back(
                {
                    id,pm(rng),x
                }

                );
            }

            auto ex=exhaustive_rate1_top_l(p,4);
            auto a=portfolio_rate1_top_l(alg,p,4);
            used.insert(a.ir.engine);
            check_candidates_equal(a.candidates,ex.candidates);
        }

        CHECK_TRUE(!used.empty());
    }

    );
}
