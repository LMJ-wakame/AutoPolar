#include "test_util.hpp"
#include <random>
using namespace polar_scl;
int main()
{
    return run_test_main("rate1",[]
    {
        {
            std::vector<ParentPath> p{{0,1.0,{0.2,-2.0,0.7}}};
            auto r=exhaustive_rate1_top_l(p,4);
            CHECK_EQ(r.candidates.size(),4u);
            CHECK_TRUE((local_rank_key(r.candidates[0])==LocalRankKey{1.0,0,0}));
            CHECK_TRUE((local_rank_key(r.candidates[1])==LocalRankKey{1.2,0,1}));
            CHECK_TRUE((local_rank_key(r.candidates[2])==LocalRankKey{1.7,0,4}));
            CHECK_TRUE((local_rank_key(r.candidates[3])==LocalRankKey{1.9,0,5}));
            CHECK_EQ(r.candidates[0].bits,(std::vector<Bit>{0,1,0}));
        }

        std::mt19937_64 rng(7);
        std::uniform_int_distribution<int> iv(-3,3), mv(0,3);
        for(std::size_t ns=1;ns<=8;++ns) for(std::size_t L:
        {
            1u,2u,4u,8u
        }

        ) for(int t=0;t<20;++t)
        {
            std::vector<ParentPath> parents;
            for(std::size_t id=0;id<L;++id)
            {
                std::vector<double> llrs(ns);
                for(double&x:llrs)x=iv(rng);
                parents.push_back(
                {
                    id,static_cast<double>(mv(rng)),llrs
                }

                );
            }

            auto ex=exhaustive_rate1_top_l(parents,L);
            auto lazy=lazy_exact_rate1_top_l(parents,L);
            auto fast=fast_sscl_bounded_rate1_top_l(parents,L);
            check_candidates_equal(lazy.candidates,ex.candidates);
            check_candidates_equal(fast.candidates,ex.candidates);
        }

        std::uniform_real_distribution<double> ur(-4,4), mr(0,2);
        std::vector<Rate1Policy> policies{{0,1,ExpansionMode::MinimumBound,8},{1,1,ExpansionMode::MinimumBound,64},{0,2,ExpansionMode::AllCompetitive,1024}};
        for(std::size_t ns=1;ns<=7;++ns) for(std::size_t L:
        {
            1u,2u,4u
        }

        ) for(int t=0;t<10;++t)
        {
            std::vector<ParentPath> p;
            for(std::size_t id=0;id<L;++id)
            {
                std::vector<double>x(ns);
                for(double&v:x)v=ur(rng);
                p.push_back(
                {
                    id,mr(rng),x
                }

                );
            }

            auto ex=exhaustive_rate1_top_l(p,L);
            for(auto pol:policies)
            {
                auto a=adaptive_rate1_top_l(p,L,pol);
                check_candidates_equal(a.candidates,ex.candidates);
            }
        }

        {
            std::vector<ParentPath> p{{0,0,{.1,.2,.3,.4}},{1,.1,{-.2,.3,-.4,.5}}};
            auto ex=exhaustive_rate1_top_l(p,4);
            auto a=adaptive_rate1_top_l(p,4,Rate1Policy{3,1,ExpansionMode::MinimumBound,2});
            CHECK_TRUE(a.stats.fallback_used);
            check_candidates_equal(a.candidates,ex.candidates);
        }

        {
            std::vector<ParentPath> p{{0,0,{.1,5,6,7,8,9}},{1,.4,{.2,5.5,6.5,7.5,8.5,9.5}}};
            auto ex=exhaustive_rate1_top_l(p,2);
            auto a=adaptive_rate1_top_l(p,2);
            check_candidates_equal(a.candidates,ex.candidates);
            CHECK_TRUE(a.stats.candidates_materialized < ex.stats.candidates_materialized);
            CHECK_TRUE(a.stats.certificate_hidden_key.has_value());
        }
    }

    );
}
