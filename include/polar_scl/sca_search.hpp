#pragma once

#include "polar_scl/algorithm_search.hpp"
#include "polar_scl/experiment.hpp"

#include <cstddef>
#include <filesystem>
#include <vector>

namespace polar_scl
{
    enum class SelectorObjective
    {
        ComplexityProxy,
        Latency
    };

    struct FrameCostWeights
    {
        double leaf_step{1.0};
        double path_split{0.25};
        double prune_event{2.0};
        double rate1_candidate{0.10};
        double frontier_expansion{0.50};
        double rate1_fallback{25.0};
        double tie_fallback{25.0};
    };

    struct SelectorScenario
    {
        Dataset dataset;
        std::vector<std::size_t> list_sizes{4,8};
    };

    struct DatasetSelectorEvaluationConfig
    {
        std::vector<SelectorScenario> scenarios;
        FrameCostWeights weights{};
        SelectorObjective objective{SelectorObjective::ComplexityProxy};
        std::size_t timing_repeats{1};
    };

    struct DatasetSelectorResult
    {
        PolarAlgorithmIR algorithm;
        std::size_t exact_cases{0};
        std::size_t total_cases{0};
        double mean_proxy_cost{0.0};
        double p99_proxy_cost{0.0};
        double mean_candidates{0.0};
        double fallback_rate{0.0};
        double mean_latency_us{0.0};
        double p99_latency_us{0.0};
        double objective_cost{0.0};
        double objective_p99{0.0};
        std::vector<std::pair<Rate1Engine,std::size_t>> engine_counts;

        bool is_exact() const
        {
            return exact_cases==total_cases;
        }
    };

    struct SCASearchConfig
    {
        SelectorSearchSpace search_space{};
        SelectorEvolutionConfig evolution{};
        std::size_t validation_top_k{12};
    };

    struct SCASearchResult
    {
        std::vector<DatasetSelectorResult> train_evaluated;
        std::vector<DatasetSelectorResult> validation_evaluated;
        DatasetSelectorResult train_best;
        DatasetSelectorResult validation_best;
        std::vector<double> best_train_cost_by_generation;
    };

    DatasetSelectorResult evaluate_selector_on_datasets(
        const SelectorGenome& genome,
        const DatasetSelectorEvaluationConfig& config);

    SCASearchResult search_sca_selector(
        const DatasetSelectorEvaluationConfig& train,
        const DatasetSelectorEvaluationConfig& validation,
        const SCASearchConfig& config={});

    void save_selector_genome(
        const SelectorGenome& genome,
        const std::filesystem::path& path);

    SelectorGenome load_selector_genome(const std::filesystem::path& path);

    void write_selector_search_csv(
        const std::vector<DatasetSelectorResult>& results,
        const std::filesystem::path& path);
}
