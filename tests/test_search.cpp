#include "test_util.hpp"
using namespace polar_scl;
int main()
{
    return run_test_main("search",[]
    {
        EvaluationConfig eval;
        eval.node_sizes={8};
        eval.list_sizes={2,4};
        eval.trials_per_case=3;
        eval.seed=5;
        EvolutionConfig evo;
        evo.population_size=4;
        evo.offspring_per_generation=3;
        evo.generations=2;
        evo.seed=9;
        auto a=evolutionary_search({},evo,eval);
        auto b=evolutionary_search({},evo,eval);
        CHECK_TRUE(!a.evaluated.empty());
        CHECK_TRUE(!a.pareto_front.empty());
        CHECK_EQ(a.evaluated.size(),b.evaluated.size());
        CHECK_EQ(a.pareto_front.size(),b.pareto_front.size());
        CHECK_EQ(a.best_mean_by_generation,b.best_mean_by_generation);
        for(const auto&r:a.evaluated)CHECK_TRUE(r.is_exact());
        SearchSpace tiny;
        tiny.initial_active_bits={0};
        tiny.growth_steps={1};
        tiny.expansion_modes={ExpansionMode::MinimumBound};
        tiny.candidate_budgets={32};
        EvolutionConfig bad;
        bad.population_size=2;
        CHECK_THROW(evolutionary_search(tiny,bad,eval));
    }

    );
}
