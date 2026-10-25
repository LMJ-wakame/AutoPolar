#include "test_util.hpp"
#include <random>
using namespace polar_scl;
int main()
{
    return run_test_main("decoder",[]
    {
        std::mt19937_64 rng(3);
        std::uniform_int_distribution<int> bit(0,1);
        for(std::size_t n:{1u,2u,4u,8u,16u,32u}) for(int t=0;t<20;++t)
        {
            std::vector<Bit>b(n);
            for(auto&x:b)x=static_cast<Bit>(bit(rng));
            CHECK_EQ(polar_transform(polar_transform(b)),b);
        }

        {
            std::vector<bool> f{true,true,false,false,true,false,false,false};
            auto ops=rate1_schedule(f);
            CHECK_EQ(ops.size(),5u);
            CHECK_TRUE(ops[0].kind==DecodeOp::Kind::Rate0&&ops[0].start==0&&ops[0].size==2);
            CHECK_TRUE(ops[1].kind==DecodeOp::Kind::Rate1&&ops[1].start==2&&ops[1].size==2);
            CHECK_TRUE(ops[4].kind==DecodeOp::Kind::Rate1&&ops[4].start==6&&ops[4].size==2);
        }

        std::normal_distribution<double> normal(0,3);
        std::bernoulli_distribution frozen_d(.45);
        std::vector<Rate1Policy> policies{{0,1,ExpansionMode::MinimumBound,64},{1,1,ExpansionMode::MinimumBound,256},{0,2,ExpansionMode::AllCompetitive,256}};
        for(std::size_t n:{4u,8u,16u}) for(std::size_t L:
        {
            1u,2u,4u,8u
        }

        ) for(int t=0;t<20;++t)
        {
            std::vector<double> llrs(n);
            for(double&v:llrs)v=normal(rng);
            std::vector<bool> f(n);
            for(std::size_t i=0;i<n;++i)f[i]=frozen_d(rng);
            auto ref=scl_decode(llrs,f,L);
            for(auto p:policies)
            {
                auto fast=rate1_scl_decode(llrs,f,L,p);
                check_paths_equal(fast.paths,ref.paths,1e-11);
            }
        }

        {
            std::vector<double> llrs{2,-1,.5,3,-2,1.5,-.7,4};
            std::vector<bool> f{true,true,false,false,false,false,false,false};
            auto ref=scl_decode(llrs,f,4);
            auto fast=rate1_scl_decode(llrs,f,4);
            check_paths_equal(fast.paths,ref.paths);
            CHECK_TRUE(fast.stats.leaf_steps<ref.stats.leaf_steps);
            CHECK_TRUE(fast.stats.rate1_nodes>0);
        }

        {
            auto f=pw_frozen_pattern(64,40);
            Rate1Policy p
            {
                0,1,ExpansionMode::MinimumBound,256
            };
            std::size_t r0=0,r1=0;
            for(int t=0;t<6;++t)
            {
                std::vector<double> llrs(64);
                for(double&v:llrs)v=normal(rng);
                auto ref=scl_decode(llrs,f,8);
                auto fast=rate1_scl_decode(llrs,f,8,p);
                check_paths_equal(fast.paths,ref.paths,1e-10);
                r0+=fast.stats.rate0_nodes;
                r1+=fast.stats.rate1_nodes;
            }

            CHECK_TRUE(r0>0&&r1>0);
        }

        std::uniform_int_distribution<int> qi(-8,8);
        Rate1Policy p
        {
            0,1,ExpansionMode::MinimumBound,256
        };
        for(std::size_t n:{4u,8u,16u})for(std::size_t L:
        {
            1u,2u,4u
        }

        )for(int t=0;t<15;++t)
        {
            std::vector<double>llrs(n);
            for(double&v:llrs)v=qi(rng);
            std::vector<bool>f(n);
            for(std::size_t i=0;i<n;++i)f[i]=frozen_d(rng);
            auto ref=scl_decode(llrs,f,L);
            auto fast=rate1_scl_decode(llrs,f,L,p);
            check_paths_equal(fast.paths,ref.paths,0.0);
        }
    }

    );
}
