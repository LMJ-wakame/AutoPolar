#include "polar_scl/search.hpp"

#include <algorithm>
#include <map>
#include <random>
#include <stdexcept>
#include <tuple>

namespace polar_scl
{
    namespace
    {
        template <class T>
        const T& choose(const std::vector<T>& values, std::mt19937_64& rng)
        {
            std::uniform_int_distribution<std::size_t> dist(0, values.size() - 1);
            return values[dist(rng)];
        }

        Rate1Policy random_policy(const SearchSpace& space, std::mt19937_64& rng)
        {
            return {choose(space.initial_active_bits, rng), choose(space.growth_steps, rng), choose(space.expansion_modes, rng), choose(space.candidate_budgets, rng)};
        }

        Rate1Policy mutate_policy(const Rate1Policy& policy, const SearchSpace& space, std::mt19937_64& rng)
        {
            std::vector<int> slots;
            if(space.initial_active_bits.size() > 1) slots.push_back(0);
            if(space.growth_steps.size() > 1) slots.push_back(1);
            if(space.expansion_modes.size() > 1) slots.push_back(2);
            if(space.candidate_budgets.size() > 1) slots.push_back(3);
            if(slots.empty()) return policy;
            std::uniform_int_distribution<std::size_t> slot_dist(0, slots.size() - 1);
            const int slot = slots[slot_dist(rng)];
            auto child = policy;
            auto choose_different = [&](const auto& values, const auto& current)
            {
                using T = std::decay_t<decltype(current)>;
                std::vector<T> choices;
                for(const auto& v : values) if(v != current) choices.push_back(v);
                return choose(choices, rng);
            };
            if(slot == 0) child.initial_active_bits = choose_different(space.initial_active_bits, policy.initial_active_bits);
            if(slot == 1) child.growth_step = choose_different(space.growth_steps, policy.growth_step);
            if(slot == 2) child.expansion_mode = choose_different(space.expansion_modes, policy.expansion_mode);
            if(slot == 3) child.candidate_budget = choose_different(space.candidate_budgets, policy.candidate_budget);
            validate_policy(child);
            return child;
        }

        auto selection_key(const PolicyResult& r)
        {
            return std::tuple{!r.is_exact(), r.fallback_rate > 0.0, r.mean_materialized, r.p99_materialized, r.mean_frontier_expansions};
        }
    }

    // namespace

    std::vector<Rate1Policy> SearchSpace::enumerate_policies() const
    {
        std::vector<Rate1Policy> result;
        for(auto initial : initial_active_bits) for(auto growth : growth_steps) for(auto mode : expansion_modes) for(auto budget : candidate_budgets)
        result.push_back({initial, growth, mode, budget});
        return result;
    }

    void SearchSpace::validate() const
    {
        if(initial_active_bits.empty() || growth_steps.empty() || expansion_modes.empty() || candidate_budgets.empty())
        throw std::invalid_argument("every search-space field must be non-empty");
        for(const auto& policy : enumerate_policies()) validate_policy(policy);
    }

    SearchResult evolutionary_search(const SearchSpace& space, const EvolutionConfig& evo, const EvaluationConfig& evaluation)
    {
        space.validate();
        auto all = space.enumerate_policies();
        if(evo.population_size < 1 || evo.offspring_per_generation < 1 || evo.generations < 1)
        throw std::invalid_argument("population/offspring/generations must be positive");
        if(evo.population_size > all.size() || evo.offspring_per_generation > all.size())
        throw std::invalid_argument("search size exceeds distinct policy space");
        std::mt19937_64 rng(evo.seed);
        std::map<Rate1Policy, PolicyResult> cache;
        auto score = [&](const Rate1Policy& p) -> const PolicyResult&
        {
            auto [it, inserted] = cache.try_emplace(p);
            if(inserted) it->second = evaluate_policy(p, evaluation);
            return it->second;
        };
        std::vector<Rate1Policy> population;
        while(population.size() < evo.population_size)
        {
            auto p = random_policy(space, rng);
            if(std::find(population.begin(), population.end(), p) == population.end()) population.push_back(p);
        }

        std::vector<double> history;
        for(std::size_t gen = 0; gen < evo.generations; ++gen)
        {
            std::vector<PolicyResult> parents;
            for(const auto& p : population) parents.push_back(score(p));
            std::sort(parents.begin(), parents.end(), [](const auto& a, const auto& b)
            {
                return selection_key(a) < selection_key(b);
            }

            );
            std::vector<Rate1Policy> children;
            std::size_t attempts = 0;
            while(children.size() < evo.offspring_per_generation && attempts++ < 20 * evo.offspring_per_generation)
            {
                std::vector<std::size_t> idx(parents.size());
                std::iota(idx.begin(), idx.end(), 0);
                std::shuffle(idx.begin(), idx.end(), rng);
                idx.resize(std::min<std::size_t>(3, idx.size()));
                auto best = idx.front();
                for(auto i : idx) if(selection_key(parents[i]) < selection_key(parents[best])) best = i;
                auto child = mutate_policy(parents[best].policy, space, rng);
                if(std::find(children.begin(), children.end(), child) == children.end()) children.push_back(child);
            }

            if(children.size() < evo.offspring_per_generation)
            {
                std::shuffle(all.begin(), all.end(), rng);
                for(const auto& p : all)
                {
                    if(children.size() >= evo.offspring_per_generation) break;
                    if(std::find(children.begin(), children.end(), p) == children.end()) children.push_back(p);
                }
            }

            std::map<Rate1Policy, PolicyResult> pool;
            for(const auto& r : parents) pool[r.policy] = r;
            for(const auto& c : children) pool[c] = score(c);
            std::vector<PolicyResult> ranked;
            for(const auto& [_, r] : pool) ranked.push_back(r);
            std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b)
            {
                return selection_key(a) < selection_key(b);
            }

            );
            population.clear();
            for(std::size_t i = 0; i < std::min(evo.population_size, ranked.size()); ++i) population.push_back(ranked[i].policy);
            double best = std::numeric_limits<double>::infinity();
            for(const auto& r : ranked) if(r.is_exact()) best = std::min(best, r.mean_materialized);
            history.push_back(best);
        }

        std::vector<PolicyResult> evaluated;
        for(const auto& [_, r] : cache) evaluated.push_back(r);
        std::sort(evaluated.begin(), evaluated.end(), [](const auto& a, const auto& b)
        {
            return selection_key(a) < selection_key(b);
        }

        );
        return {evaluated, pareto_policies(evaluated), history};
    }
}

// namespace polar_scl
