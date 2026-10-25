#include "polar_scl/construction.hpp"
#include "polar_scl/portfolio.hpp"
#include "polar_scl/special_nodes.hpp"
#include "test_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

using namespace polar_scl;

namespace
{
    std::vector<LocalRankKey> keys(const std::vector<Rate1Candidate>& candidates)
    {
        std::vector<LocalRankKey> result;
        for(const auto& candidate:candidates)
        {
            result.push_back(local_rank_key(candidate));
        }
        return result;
    }

    Rate1SearchResult exhaustive_spc(const std::vector<ParentPath>& parents,std::size_t list_size)
    {
        std::vector<Rate1Candidate> candidates;
        for(const auto& parent:parents)
        {
            const std::size_t n=parent.llrs.size();
            const std::uint64_t total=1ULL<<n;
            for(std::uint64_t mask=0;mask<total;++mask)
            {
                std::vector<Bit> bits(n);
                double penalty=0.0;
                unsigned parity=0;
                for(std::size_t i=0;i<n;++i)
                {
                    Bit bit=parent.llrs[i]<0.0?Bit{1}:Bit{0};
                    if(((mask>>i)&1ULL)!=0ULL)
                    {
                        bit^=Bit{1};
                        penalty+=std::abs(parent.llrs[i]);
                    }
                    bits[i]=bit;
                    parity^=static_cast<unsigned>(bit);
                }
                if(parity==0)
                {
                    candidates.push_back({parent.parent_id,parent.path_metric+penalty,penalty,mask,std::move(bits)});
                }
            }
        }
        std::sort(candidates.begin(),candidates.end(),[](const auto& lhs,const auto& rhs)
        {
            return local_rank_key(lhs)<local_rank_key(rhs);
        });
        if(candidates.size()>list_size)
        {
            candidates.resize(list_size);
        }
        return {std::move(candidates),{}};
    }
}

int main()
{
    std::mt19937_64 rng(20260816);
    std::normal_distribution<double> llr(0.0,2.0);
    for(std::size_t node_size:{2U,4U,8U})
    {
        for(std::size_t list_size:{2U,4U,8U})
        {
            for(int trial=0;trial<80;++trial)
            {
                std::vector<ParentPath> parents;
                const std::size_t parent_count=std::min<std::size_t>(list_size,4);
                for(std::size_t p=0;p<parent_count;++p)
                {
                    std::vector<double> values(node_size);
                    for(double& value:values)
                    {
                        value=std::round(llr(rng)*4.0)/4.0;
                    }
                    parents.push_back({p,0.1*static_cast<double>(p),std::move(values)});
                }
                const auto expected=exhaustive_spc(parents,list_size);
                const auto actual=spc_node_top_l(parents,list_size);
                CHECK_TRUE(keys(expected.candidates)==keys(actual.candidates));

                const auto rep=repetition_node_top_l(parents,list_size);
                CHECK_TRUE(!rep.candidates.empty());
                for(const auto& candidate:rep.candidates)
                {
                    CHECK_TRUE(std::all_of(candidate.bits.begin(),candidate.bits.end(),[&](Bit bit)
                    {
                        return bit==candidate.bits.front();
                    }));
                }
            }
        }
    }

    for(std::size_t n:{64U,128U})
    {
        const std::size_t k=n/2;
        const auto frozen=pw_frozen_pattern(n,k);
        SelectorGenome genome;
        genome.enable_repetition=true;
        genome.enable_spc=true;
        genome.list_size_threshold=4;
        genome.low_list_engine=Rate1Engine::CertifiedFrontier;
        genome.small_engine=Rate1Engine::FastSSCLBounded;
        genome.easy_engine=Rate1Engine::CertifiedFrontier;
        genome.hard_engine=Rate1Engine::FastSSCLBounded;
        PolarAlgorithmIR algorithm;
        algorithm.selector_mode=SelectorMode::Rule;
        algorithm.selector=genome;
        std::normal_distribution<double> channel(0.0,3.0);
        for(int trial=0;trial<80;++trial)
        {
            std::vector<double> llrs(n);
            for(double& value:llrs)
            {
                value=std::round(channel(rng)*8.0)/8.0;
            }
            const auto reference=scl_decode(llrs,frozen,8);
            const auto actual=portfolio_scl_decode(llrs,frozen,8,algorithm);
            CHECK_TRUE(paths_equivalent(reference.paths,actual.paths));
        }
    }
    return 0;
}
