#include "polar_scl/special_nodes.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <tuple>

namespace polar_scl
{
    namespace
    {
        GeneralRankKey rank_candidate(const Rate1Candidate& candidate,const CandidateRankFn& rank_key)
        {
            return rank_key?rank_key(candidate):general_local_rank_key(candidate);
        }

        bool rank_less(const GeneralRankKey& lhs,const GeneralRankKey& rhs)
        {
            if(lhs.metric<rhs.metric)
            {
                return true;
            }
            if(rhs.metric<lhs.metric)
            {
                return false;
            }
            return std::lexicographical_compare(lhs.tie.begin(),lhs.tie.end(),rhs.tie.begin(),rhs.tie.end());
        }

        std::vector<Bit> hard_bits(const std::vector<double>& llrs)
        {
            std::vector<Bit> bits;
            bits.reserve(llrs.size());
            for(double llr:llrs)
            {
                bits.push_back(llr<0.0?Bit{1}:Bit{0});
            }
            return bits;
        }

        Rate1Candidate make_flip_candidate(const ParentPath& parent,std::uint64_t flip_mask)
        {
            if(parent.llrs.size()>64)
            {
                throw std::invalid_argument("special-node mask implementation supports node sizes up to 64");
            }
            auto bits=hard_bits(parent.llrs);
            double penalty=0.0;
            for(std::size_t index=0;index<parent.llrs.size();++index)
            {
                if(((flip_mask>>index)&1ULL)!=0ULL)
                {
                    bits[index]^=Bit{1};
                    penalty+=std::abs(parent.llrs[index]);
                }
            }
            return {parent.parent_id,parent.path_metric+penalty,penalty,flip_mask,std::move(bits)};
        }

        class SubsetStream
        {
            public:
            SubsetStream(const ParentPath& parent,SearchStats& stats)
            :parent_(parent),stats_(stats)
            {
                positions_.resize(parent.llrs.size());
                std::iota(positions_.begin(),positions_.end(),std::size_t{0});
                std::sort(positions_.begin(),positions_.end(),[&](std::size_t lhs,std::size_t rhs)
                {
                    const double a=std::abs(parent.llrs[lhs]);
                    const double b=std::abs(parent.llrs[rhs]);
                    return a<b||(a==b&&lhs<rhs);
                });
                weights_.reserve(positions_.size());
                for(std::size_t position:positions_)
                {
                    weights_.push_back(std::abs(parent.llrs[position]));
                }
            }

            std::optional<Rate1Candidate> next_matching(unsigned parity)
            {
                while(true)
                {
                    auto candidate=next_any();
                    if(!candidate)
                    {
                        return std::nullopt;
                    }
                    if((std::popcount(candidate->flip_mask)&1U)==parity)
                    {
                        return candidate;
                    }
                }
            }

            private:
            struct State
            {
                double penalty{};
                std::uint64_t original_mask{};
                std::size_t last_index{};
                std::uint64_t sorted_mask{};
            };

            struct Greater
            {
                bool operator()(const State& lhs,const State& rhs) const
                {
                    return std::tie(lhs.penalty,lhs.original_mask,lhs.last_index,lhs.sorted_mask)>
                        std::tie(rhs.penalty,rhs.original_mask,rhs.last_index,rhs.sorted_mask);
                }
            };

            std::optional<Rate1Candidate> next_any()
            {
                if(!yielded_empty_)
                {
                    yielded_empty_=true;
                    return make_flip_candidate(parent_,0);
                }
                if(weights_.empty())
                {
                    return std::nullopt;
                }
                if(!initialized_)
                {
                    initialized_=true;
                    push(weights_[0],1ULL,0);
                }
                if(heap_.empty())
                {
                    return std::nullopt;
                }
                const State state=heap_.top();
                heap_.pop();
                ++stats_.subset_states_popped;
                const std::size_t next_index=state.last_index+1;
                if(next_index<weights_.size())
                {
                    const std::uint64_t next_bit=1ULL<<next_index;
                    push(state.penalty+weights_[next_index],state.sorted_mask|next_bit,next_index);
                    push(state.penalty-weights_[state.last_index]+weights_[next_index],
                        (state.sorted_mask^(1ULL<<state.last_index))|next_bit,next_index);
                }
                return make_flip_candidate(parent_,state.original_mask);
            }

            std::uint64_t to_original_mask(std::uint64_t sorted_mask) const
            {
                std::uint64_t result=0;
                for(std::size_t bit=0;bit<positions_.size();++bit)
                {
                    if(((sorted_mask>>bit)&1ULL)!=0ULL)
                    {
                        result|=1ULL<<positions_[bit];
                    }
                }
                return result;
            }

            void push(double penalty,std::uint64_t sorted_mask,std::size_t last_index)
            {
                heap_.push({penalty,to_original_mask(sorted_mask),last_index,sorted_mask});
                ++stats_.subset_states_pushed;
            }

            const ParentPath& parent_;
            SearchStats& stats_;
            std::vector<std::size_t> positions_;
            std::vector<double> weights_;
            bool yielded_empty_{false};
            bool initialized_{false};
            std::priority_queue<State,std::vector<State>,Greater> heap_;
        };
    }

    Rate1SearchResult repetition_node_top_l(
        const std::vector<ParentPath>& parents,
        std::size_t list_size,
        CandidateRankFn rank_key)
    {
        validate_parent_paths(parents,list_size);
        SearchStats stats;
        stats.algorithm="repetition_exact";
        std::vector<Rate1Candidate> candidates;
        candidates.reserve(parents.size()*2);
        for(const ParentPath& parent:parents)
        {
            if(parent.llrs.size()>64)
            {
                throw std::invalid_argument("repetition node size above 64 is not supported");
            }
            for(Bit value:{Bit{0},Bit{1}})
            {
                std::vector<Bit> bits(parent.llrs.size(),value);
                double penalty=0.0;
                std::uint64_t mask=0;
                for(std::size_t index=0;index<parent.llrs.size();++index)
                {
                    const Bit hard=parent.llrs[index]<0.0?Bit{1}:Bit{0};
                    if(hard!=value)
                    {
                        penalty+=std::abs(parent.llrs[index]);
                        mask|=1ULL<<index;
                    }
                }
                candidates.push_back({parent.parent_id,parent.path_metric+penalty,penalty,mask,std::move(bits)});
            }
        }
        stats.candidates_materialized=candidates.size();
        std::sort(candidates.begin(),candidates.end(),[&](const Rate1Candidate& lhs,const Rate1Candidate& rhs)
        {
            return rank_less(rank_candidate(lhs,rank_key),rank_candidate(rhs,rank_key));
        });
        if(candidates.size()>list_size)
        {
            candidates.resize(list_size);
        }
        return {std::move(candidates),std::move(stats)};
    }

    Rate1SearchResult spc_node_top_l(
        const std::vector<ParentPath>& parents,
        std::size_t list_size,
        CandidateRankFn rank_key)
    {
        validate_parent_paths(parents,list_size);
        if(parents.front().llrs.size()<2)
        {
            throw std::invalid_argument("SPC node size must be at least two");
        }
        SearchStats stats;
        stats.algorithm="spc_lazy_exact";
        std::vector<SubsetStream> streams;
        streams.reserve(parents.size());
        for(const ParentPath& parent:parents)
        {
            streams.emplace_back(parent,stats);
        }

        struct HeapItem
        {
            Rate1Candidate candidate;
            GeneralRankKey key;
            std::size_t stream_index{};
        };
        struct Greater
        {
            bool operator()(const HeapItem& lhs,const HeapItem& rhs) const
            {
                if(lhs.key.metric!=rhs.key.metric)
                {
                    return lhs.key.metric>rhs.key.metric;
                }
                return std::lexicographical_compare(
                    rhs.key.tie.begin(),rhs.key.tie.end(),lhs.key.tie.begin(),lhs.key.tie.end());
            }
        };

        std::priority_queue<HeapItem,std::vector<HeapItem>,Greater> heap;
        std::vector<unsigned> required_parity(parents.size(),0);
        for(std::size_t index=0;index<parents.size();++index)
        {
            unsigned hard_parity=0;
            for(Bit bit:hard_bits(parents[index].llrs))
            {
                hard_parity^=static_cast<unsigned>(bit);
            }
            required_parity[index]=hard_parity;
            auto candidate=streams[index].next_matching(hard_parity);
            if(candidate)
            {
                heap.push({*candidate,rank_candidate(*candidate,rank_key),index});
                ++stats.global_heap_pushes;
            }
        }

        std::vector<Rate1Candidate> output;
        output.reserve(list_size);
        while(!heap.empty()&&output.size()<list_size)
        {
            HeapItem item=heap.top();
            heap.pop();
            ++stats.global_heap_pops;
            output.push_back(std::move(item.candidate));
            const std::size_t index=item.stream_index;
            auto next=streams[index].next_matching(required_parity[index]);
            if(next)
            {
                heap.push({*next,rank_candidate(*next,rank_key),index});
                ++stats.global_heap_pushes;
            }
        }
        stats.candidates_materialized=output.size();
        return {std::move(output),std::move(stats)};
    }
}
