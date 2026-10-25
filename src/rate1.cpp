#include "polar_scl/rate1.hpp"
#include "polar_scl_lab/linear.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <utility>

namespace polar_scl
{
    namespace
    {
        bool rank_less(const GeneralRankKey& lhs, const GeneralRankKey& rhs)
        {
            if(lhs.metric < rhs.metric) return true;
            if(rhs.metric < lhs.metric) return false;
            return std::lexicographical_compare(lhs.tie.begin(), lhs.tie.end(), rhs.tie.begin(), rhs.tie.end());
        }

        Rate1Candidate make_candidate(const ParentPath& parent, std::uint64_t flip_mask)
        {
            if(parent.llrs.size() > 64)
            {
                throw std::invalid_argument("Rate-1 mask implementation supports node sizes up to 64");
            }

            return polar_scl_lab::make_candidate(parent, flip_mask);
        }

        std::size_t saturated_pow2(std::size_t exponent)
        {
            if(exponent >= std::numeric_limits<std::size_t>::digits)
            {
                return std::numeric_limits<std::size_t>::max();
            }

            return std::size_t{1} << exponent;
        }

        void merge_stats(SearchStats& target, const SearchStats& fallback)
        {
            target.subset_states_pushed += fallback.subset_states_pushed;
            target.subset_states_popped += fallback.subset_states_popped;
            target.global_heap_pushes += fallback.global_heap_pushes;
            target.global_heap_pops += fallback.global_heap_pops;
            target.candidates_materialized += fallback.candidates_materialized;
        }

        class SubsetStream
        {
            public:
            SubsetStream(const ParentPath& parent, SearchStats& stats)
            : parent_(parent), stats_(stats)
            {
                positions_.resize(parent.llrs.size());
                std::iota(positions_.begin(), positions_.end(), std::size_t{0}

                );
                std::sort(positions_.begin(), positions_.end(), [&](std::size_t a, std::size_t b)
                {
                    const auto aa = std::abs(parent.llrs[a]);
                    const auto bb = std::abs(parent.llrs[b]);
                    return aa < bb || (aa == bb && a < b);
                }

                );
                weights_.reserve(positions_.size());
                for(auto pos : positions_) weights_.push_back(std::abs(parent.llrs[pos]));
            }

            std::optional<Rate1Candidate> next_candidate()
            {
                if(!yielded_empty_)
                {
                    yielded_empty_ = true;
                    return make_candidate(parent_, 0);
                }

                if(!initialized_)
                {
                    initialized_ = true;
                    push(weights_.at(0), 1ULL, 0);
                }

                if(heap_.empty()) return std::nullopt;

                auto state = heap_.top();
                heap_.pop();
                ++stats_.subset_states_popped;

                const std::size_t next_index = state.last_index + 1;
                if(next_index < weights_.size())
                {
                    const std::uint64_t next_bit = 1ULL << next_index;
                    push(state.penalty + weights_[next_index], state.sorted_mask | next_bit, next_index);
                    push(state.penalty - weights_[state.last_index] + weights_[next_index],
                    (state.sorted_mask ^ (1ULL << state.last_index)) | next_bit,
                    next_index);
                }

                return make_candidate(parent_, state.original_mask);
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
                bool operator()(const State& a, const State& b) const
                {
                    return std::tie(a.penalty, a.original_mask, a.last_index, a.sorted_mask) >
                    std::tie(b.penalty, b.original_mask, b.last_index, b.sorted_mask);
                }
            };

            std::uint64_t original_mask(std::uint64_t sorted_mask) const
            {
                std::uint64_t result = 0;
                for(std::size_t bit = 0; bit < positions_.size(); ++bit)
                {
                    if(((sorted_mask >> bit) & 1ULL) != 0ULL) result |= 1ULL << positions_[bit];
                }

                return result;
            }

            void push(double penalty, std::uint64_t sorted_mask, std::size_t last_index)
            {
                heap_.push(State{penalty, original_mask(sorted_mask), last_index, sorted_mask});
                ++stats_.subset_states_pushed;
            }

            const ParentPath& parent_;
            SearchStats& stats_;
            std::vector<std::size_t> positions_;
            std::vector<double> weights_;
            bool yielded_empty_{false};
            bool initialized_{false};
            std::priority_queue<State, std::vector<State>, Greater> heap_;
        };

        class BoundedSubsetStream
        {
            public:
            BoundedSubsetStream(const ParentPath& parent,std::size_t active_bits,SearchStats& stats)
            : parent_(parent),stats_(stats)
            {
                positions_.resize(parent.llrs.size());
                std::iota(positions_.begin(),positions_.end(),std::size_t{0}

                );
                std::sort(positions_.begin(),positions_.end(),[&](std::size_t a,std::size_t b)
                {
                    const auto aa=std::abs(parent.llrs[a]);
                    const auto bb=std::abs(parent.llrs[b]);
                    return aa<bb||(aa==bb&&a<b);
                }
                );
                positions_.resize(std::min(active_bits,positions_.size()));
                weights_.reserve(positions_.size());
                for(auto position:positions_)
                {
                    weights_.push_back(std::abs(parent.llrs[position]));
                }
            }

            std::optional<Rate1Candidate> next_candidate()
            {
                if(!yielded_empty_)
                {
                    yielded_empty_=true;
                    return make_candidate(parent_,0);
                }

                if(weights_.empty())
                {
                    return std::nullopt;
                }

                if(!initialized_)
                {
                    initialized_=true;
                    push(weights_.front(),1ULL,0);
                }

                if(heap_.empty())
                {
                    return std::nullopt;
                }

                auto state=heap_.top();
                heap_.pop();
                ++stats_.subset_states_popped;

                const std::size_t next_index=state.last_index+1;
                if(next_index<weights_.size())
                {
                    const std::uint64_t next_bit=1ULL<<next_index;
                    push(state.penalty+weights_[next_index],state.sorted_mask|next_bit,next_index);
                    push(
                    state.penalty-weights_[state.last_index]+weights_[next_index],
                    (state.sorted_mask^(1ULL<<state.last_index))|next_bit,
                    next_index);
                }

                return make_candidate(parent_,state.original_mask);
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
                bool operator()(const State& a,const State& b) const
                {
                    return std::tie(a.penalty,a.original_mask,a.last_index,a.sorted_mask)>
                    std::tie(b.penalty,b.original_mask,b.last_index,b.sorted_mask);
                }
            };

            std::uint64_t original_mask(std::uint64_t sorted_mask) const
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
                heap_.push(State{penalty,original_mask(sorted_mask),last_index,sorted_mask});
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

        struct FrontierState
        {
            const ParentPath* parent{};
            std::vector<std::size_t> positions;
            std::vector<double> weights;
            std::size_t active_bits{};

            static FrontierState create(const ParentPath& parent, std::size_t active)
            {
                FrontierState state;
                state.parent = &parent;
                state.positions.resize(parent.llrs.size());
                std::iota(state.positions.begin(), state.positions.end(), std::size_t{0}

                );
                std::sort(state.positions.begin(), state.positions.end(), [&](std::size_t a, std::size_t b)
                {
                    const auto aa = std::abs(parent.llrs[a]);
                    const auto bb = std::abs(parent.llrs[b]);
                    return aa < bb || (aa == bb && a < b);
                }

                );
                for(auto pos : state.positions) state.weights.push_back(std::abs(parent.llrs[pos]));
                state.active_bits = std::min(active, state.positions.size());
                return state;
            }

            std::uint64_t original_mask(std::uint64_t sorted_mask) const
            {
                std::uint64_t mask = 0;
                for(std::size_t i = 0; i < active_bits; ++i)
                {
                    if(((sorted_mask >> i) & 1ULL) != 0ULL) mask |= 1ULL << positions[i];
                }

                return mask;
            }

            std::optional<double> hidden_lower_metric() const
            {
                if(active_bits >= positions.size()) return std::nullopt;
                return parent->path_metric + weights[active_bits];
            }

            std::optional<GeneralRankKey> hidden_lower_key() const
            {
                const auto metric = hidden_lower_metric();
                if(!metric) return std::nullopt;
                const double hidden_weight = weights[active_bits];
                std::size_t smallest_position = positions[active_bits];
                for(std::size_t i = active_bits; i < positions.size(); ++i)
                {
                    if(weights[i] == hidden_weight) smallest_position = std::min(smallest_position, positions[i]);
                }

                return GeneralRankKey{*metric, {static_cast<std::uint64_t>(parent->parent_id), 1ULL << smallest_position}};
            }
        };

        std::vector<Rate1Candidate> frontier_candidates(const FrontierState& state, std::optional<std::size_t> old_active)
        {
            const std::size_t start = old_active ? saturated_pow2(*old_active) : 0;
            const std::size_t stop = saturated_pow2(state.active_bits);
            if(stop == std::numeric_limits<std::size_t>::max())
            {
                throw std::overflow_error("frontier too large to enumerate");
            }

            std::vector<Rate1Candidate> result;
            result.reserve(stop - start);
            for(std::size_t mask = start; mask < stop; ++mask)
            {
                result.push_back(make_candidate(*state.parent, state.original_mask(static_cast<std::uint64_t>(mask))));
            }

            return result;
        }
    }

    // namespace

    bool operator<(const GeneralRankKey& lhs, const GeneralRankKey& rhs)
    {
        return rank_less(lhs, rhs);
    }

    bool operator>(const GeneralRankKey& lhs, const GeneralRankKey& rhs)
    {
        return rhs < lhs;
    }

    bool operator<=(const GeneralRankKey& lhs, const GeneralRankKey& rhs)
    {
        return !(rhs < lhs);
    }

    bool operator==(const GeneralRankKey& lhs, const GeneralRankKey& rhs)
    {
        return lhs.metric == rhs.metric && lhs.tie == rhs.tie;
    }

    LocalRankKey local_rank_key(const Rate1Candidate& candidate)
    {
        return {candidate.path_metric, candidate.parent_id, candidate.flip_mask};
    }

    GeneralRankKey general_local_rank_key(const Rate1Candidate& candidate)
    {
        return {candidate.path_metric,
            {
                static_cast<std::uint64_t>(candidate.parent_id), candidate.flip_mask
            }
        };
    }

    void validate_parent_paths(const std::vector<ParentPath>& parents, std::size_t list_size)
    {
        if(parents.empty()) throw std::invalid_argument("at least one parent path is required");
        if(list_size < 1) throw std::invalid_argument("list_size must be positive");
        const auto node_size = parents.front().llrs.size();
        if(node_size == 0 || node_size > 64) throw std::invalid_argument("Rate-1 node size must lie in [1,64]");
        std::set<std::size_t> ids;
        for(const auto& parent : parents)
        {
            if(parent.llrs.size() != node_size) throw std::invalid_argument("all parents must share node size");
            if(!std::isfinite(parent.path_metric)) throw std::invalid_argument("path metric must be finite");
            for(double value : parent.llrs)
            {
                if(!std::isfinite(value)) throw std::invalid_argument("all LLRs must be finite");
            }

            if(!ids.insert(parent.parent_id).second) throw std::invalid_argument("parent_id values must be unique");
        }
    }

    bool rate1_metric_tie_risk(const std::vector<ParentPath>& parents,double tolerance)
    {
        if(parents.empty())
        {
            return false;
        }
        auto close=[](double lhs,double rhs,double epsilon)
        {
            return std::abs(lhs-rhs)<=epsilon*std::max({1.0,std::abs(lhs),std::abs(rhs)});
        };

        std::vector<double> parent_metrics;
        parent_metrics.reserve(parents.size());
        for(const ParentPath& parent:parents)
        {
            parent_metrics.push_back(parent.path_metric);
            std::vector<double> reliabilities;
            reliabilities.reserve(parent.llrs.size());
            for(double llr:parent.llrs)
            {
                const double reliability=std::abs(llr);
                if(reliability<=tolerance)
                {
                    return true;
                }
                reliabilities.push_back(reliability);
            }
            std::sort(reliabilities.begin(),reliabilities.end());
            for(std::size_t index=1;index<reliabilities.size();++index)
            {
                if(close(reliabilities[index-1],reliabilities[index],tolerance))
                {
                    return true;
                }
            }
        }

        std::sort(parent_metrics.begin(),parent_metrics.end());
        for(std::size_t index=1;index<parent_metrics.size();++index)
        {
            if(close(parent_metrics[index-1],parent_metrics[index],tolerance))
            {
                return true;
            }
        }
        return false;
    }

    Rate1SearchResult exhaustive_rate1_top_l(const std::vector<ParentPath>& parents, std::size_t list_size)
    {
        validate_parent_paths(parents, list_size);
        const auto node_size = parents.front().llrs.size();
        if(node_size >= std::numeric_limits<std::size_t>::digits)
        {
            throw std::invalid_argument("exhaustive enumeration is unavailable for this node size");
        }

        const std::size_t per_parent = std::size_t{1} << node_size;
        std::vector<Rate1Candidate> candidates;
        if(per_parent <= std::numeric_limits<std::size_t>::max() / parents.size())
        candidates.reserve(per_parent * parents.size());
        for(const auto& parent : parents)
        {
            for(std::size_t mask = 0; mask < per_parent; ++mask)
            candidates.push_back(make_candidate(parent, static_cast<std::uint64_t>(mask)));
        }

        std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b)
        {
            return local_rank_key(a) < local_rank_key(b);
        }

        );
        SearchStats stats;
        stats.algorithm = "exhaustive";
        stats.candidates_materialized = candidates.size();
        if(candidates.size() > list_size) candidates.resize(list_size);
        return {std::move(candidates), std::move(stats)};
    }

    Rate1SearchResult lazy_exact_rate1_top_l(const std::vector<ParentPath>& parents, std::size_t list_size)
    {
        validate_parent_paths(parents, list_size);
        SearchStats stats;
        stats.algorithm = "lazy_exact";
        std::vector<SubsetStream> streams;
        streams.reserve(parents.size());
        for(const auto& parent : parents) streams.emplace_back(parent, stats);

        struct HeapItem
        {
            LocalRankKey key;
            std::size_t stream_index{};
            Rate1Candidate candidate;
        };
        struct Greater
        {
            bool operator()(const HeapItem& a, const HeapItem& b) const
            {
                return a.key > b.key;
            }
        };
        std::priority_queue<HeapItem, std::vector<HeapItem>, Greater> heap;

        for(std::size_t i = 0; i < streams.size(); ++i)
        {
            auto first = streams[i].next_candidate();
            if(first)
            {
                heap.push({local_rank_key(*first), i, *first});
                ++stats.global_heap_pushes;
                ++stats.candidates_materialized;
            }
        }

        std::vector<Rate1Candidate> selected;
        selected.reserve(list_size);
        while(!heap.empty() && selected.size() < list_size)
        {
            auto item = heap.top();
            heap.pop();
            ++stats.global_heap_pops;
            selected.push_back(std::move(item.candidate));
            auto next = streams[item.stream_index].next_candidate();
            if(next)
            {
                heap.push({local_rank_key(*next), item.stream_index, *next});
                ++stats.global_heap_pushes;
                ++stats.candidates_materialized;
            }
        }

        std::sort(selected.begin(), selected.end(), [](const auto& a, const auto& b)
        {
            return local_rank_key(a) < local_rank_key(b);
        }

        );
        return {std::move(selected), std::move(stats)};
    }

    Rate1SearchResult fast_sscl_bounded_rate1_top_l(
    const std::vector<ParentPath>& parents,
    std::size_t list_size)
    {
        validate_parent_paths(parents,list_size);
        const std::size_t node_size=parents.front().llrs.size();
        const std::size_t active=std::min(list_size-1,node_size);

        SearchStats stats;
        stats.algorithm="fast_sscl_bounded";
        std::vector<BoundedSubsetStream> streams;
        streams.reserve(parents.size());
        for(const auto& parent:parents)
        {
            streams.emplace_back(parent,active,stats);
        }

        struct HeapItem
        {
            LocalRankKey key;
            std::size_t stream_index{};
            Rate1Candidate candidate;
        };

        struct Greater
        {
            bool operator()(const HeapItem& a,const HeapItem& b) const
            {
                return a.key>b.key;
            }
        };

        std::priority_queue<HeapItem,std::vector<HeapItem>,Greater> heap;
        for(std::size_t index=0;index<streams.size();++index)
        {
            auto first=streams[index].next_candidate();
            if(first)
            {
                heap.push({local_rank_key(*first),index,*first});
                ++stats.global_heap_pushes;
                ++stats.candidates_materialized;
            }
        }

        std::vector<Rate1Candidate> selected;
        selected.reserve(list_size);
        while(!heap.empty()&&selected.size()<list_size)
        {
            auto item=heap.top();
            heap.pop();
            ++stats.global_heap_pops;
            selected.push_back(std::move(item.candidate));
            auto next=streams[item.stream_index].next_candidate();
            if(next)
            {
                heap.push({local_rank_key(*next),item.stream_index,*next});
                ++stats.global_heap_pushes;
                ++stats.candidates_materialized;
            }
        }

        std::sort(selected.begin(),selected.end(),[](const auto& a,const auto& b)
        {
            return local_rank_key(a)<local_rank_key(b);
        }
        );
        return {std::move(selected),std::move(stats)};
    }

    Rate1SearchResult adaptive_rate1_top_l(
    const std::vector<ParentPath>& parents,
    std::size_t list_size,
    const Rate1Policy& policy,
    CandidateRankFn rank_key,
    HiddenRankFn hidden_rank_key)
    {
        validate_policy(policy);
        validate_parent_paths(parents, list_size);
        if(!rank_key) rank_key = general_local_rank_key;

        SearchStats stats;
        stats.algorithm = "adaptive_certified";
        std::vector<FrontierState> states;
        states.reserve(parents.size());
        for(const auto& parent : parents) states.push_back(FrontierState::create(parent, policy.initial_active_bits));

        auto hidden_rank = [&](const FrontierState& state) -> std::optional<GeneralRankKey>
        {
            const auto metric = state.hidden_lower_metric();
            if(!metric) return std::nullopt;
            if(hidden_rank_key) return hidden_rank_key(*state.parent, *metric);
            return state.hidden_lower_key();
        };

        std::size_t initial_count = 0;
        for(const auto& state : states)
        {
            const auto count = saturated_pow2(state.active_bits);
            if(count == std::numeric_limits<std::size_t>::max() || initial_count > policy.candidate_budget - std::min(count, policy.candidate_budget))
            {
                auto fallback = lazy_exact_rate1_top_l(parents, list_size);
                stats.fallback_used = true;
                stats.fallback_reason = "initial frontier exceeds candidate budget";
                merge_stats(stats, fallback.stats);
                return {std::move(fallback.candidates), std::move(stats)};
            }

            initial_count += count;
        }

        if(initial_count > policy.candidate_budget)
        {
            auto fallback = lazy_exact_rate1_top_l(parents, list_size);
            stats.fallback_used = true;
            stats.fallback_reason = "initial frontier exceeds candidate budget";
            merge_stats(stats, fallback.stats);
            return {std::move(fallback.candidates), std::move(stats)};
        }

        std::vector<Rate1Candidate> visible;
        visible.reserve(initial_count);
        for(const auto& state : states)
        {
            auto generated = frontier_candidates(state, std::nullopt);
            visible.insert(visible.end(), generated.begin(), generated.end());
            stats.per_parent_active_bits[state.parent->parent_id] = state.active_bits;
        }

        stats.candidates_materialized = visible.size();
        const auto per_parent_total = saturated_pow2(parents.front().llrs.size());
        std::size_t total_candidates = std::numeric_limits<std::size_t>::max();
        if(per_parent_total != std::numeric_limits<std::size_t>::max() &&
        per_parent_total <= std::numeric_limits<std::size_t>::max() / parents.size())
        {
            total_candidates = per_parent_total * parents.size();
        }

        const std::size_t target = std::min(list_size, total_candidates);

        while(true)
        {
            std::sort(visible.begin(), visible.end(), [&](const auto& a, const auto& b)
            {
                return rank_key(a) < rank_key(b);
            }

            );
            const std::size_t selected_count = std::min(target, visible.size());

            std::vector<std::pair<std::size_t, GeneralRankKey>> hidden;
            for(std::size_t i = 0; i < states.size(); ++i)
            {
                if(auto key = hidden_rank(states[i])) hidden.emplace_back(i, *key);
            }

            ++stats.certificate_checks;

            GeneralRankKey worst{std::numeric_limits<double>::infinity(),
                {
                }
            };
            if(selected_count == target)
            {
                worst = rank_key(visible[target - 1]);
                stats.certificate_worst_key = worst;
                std::optional<GeneralRankKey> best_hidden;
                for(const auto& [_, key] : hidden) if(!best_hidden || key < *best_hidden) best_hidden = key;
                stats.certificate_hidden_key = best_hidden;
                if(!best_hidden || *best_hidden > worst)
                {
                    visible.resize(target);
                    return {std::move(visible), std::move(stats)};
                }
            }

            std::vector<std::size_t> competitive;
            for(const auto& [index, key] : hidden) if(key <= worst) competitive.push_back(index);
            if(competitive.empty())
            {
                auto fallback = lazy_exact_rate1_top_l(parents, list_size);
                stats.fallback_used = true;
                stats.fallback_reason = "no expandable frontier before certificate";
                merge_stats(stats, fallback.stats);
                return {std::move(fallback.candidates), std::move(stats)};
            }

            if(policy.expansion_mode == ExpansionMode::MinimumBound)
            {
                auto best = *std::min_element(competitive.begin(), competitive.end(), [&](std::size_t a, std::size_t b)
                {
                    return *hidden_rank(states[a]) < *hidden_rank(states[b]);
                }

                );
                competitive.clear();
                competitive.push_back(best);
            }

            struct Plan
            {
                std::size_t index, old_active, new_active, count;
            };
            std::vector<Plan> plans;
            std::size_t additional = 0;
            bool overflow = false;
            for(auto index : competitive)
            {
                auto& state = states[index];
                const auto old_active = state.active_bits;
                const auto new_active = std::min(old_active + policy.growth_step, state.positions.size());
                const auto old_count = saturated_pow2(old_active);
                const auto new_count = saturated_pow2(new_active);
                if(new_count == std::numeric_limits<std::size_t>::max())
                {
                    overflow = true;
                    break;
                }

                const auto count = new_count - old_count;
                if(additional > policy.candidate_budget - std::min(count, policy.candidate_budget))
                {
                    overflow = true;
                    break;
                }

                additional += count;
                plans.push_back({index, old_active, new_active, count});
            }

            if(overflow || stats.candidates_materialized + additional > policy.candidate_budget)
            {
                auto fallback = lazy_exact_rate1_top_l(parents, list_size);
                stats.fallback_used = true;
                stats.fallback_reason = "frontier expansion exceeds candidate budget";
                merge_stats(stats, fallback.stats);
                return {std::move(fallback.candidates), std::move(stats)};
            }

            for(const auto& plan : plans)
            {
                auto& state = states[plan.index];
                state.active_bits = plan.new_active;
                auto generated = frontier_candidates(state, plan.old_active);
                visible.insert(visible.end(), generated.begin(), generated.end());
                ++stats.frontier_expansions;
                stats.per_parent_active_bits[state.parent->parent_id] = state.active_bits;
            }

            stats.candidates_materialized += additional;
        }
    }
}

// namespace polar_scl
