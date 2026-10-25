#include "test_util.hpp"

#include <filesystem>

using namespace polar_scl;

int main()
{
    return run_test_main("sca_search",[]
    {
        AwgnDatasetConfig train_config;
        train_config.length=32;
        train_config.payload_bits=8;
        train_config.crc_bits=8;
        train_config.ebn0_db=1.0;
        train_config.frames=3;
        train_config.seed=101;

        AwgnDatasetConfig validation_config=train_config;
        validation_config.seed=202;

        DatasetSelectorEvaluationConfig train;
        train.scenarios.push_back({generate_awgn_dataset(train_config),{2,4}});
        DatasetSelectorEvaluationConfig validation;
        validation.scenarios.push_back({generate_awgn_dataset(validation_config),{2,4}});

        SCASearchConfig config;
        config.evolution.population_size=4;
        config.evolution.offspring_per_generation=4;
        config.evolution.generations=2;
        config.evolution.seed=303;
        config.validation_top_k=4;

        const SCASearchResult first=search_sca_selector(train,validation,config);
        const SCASearchResult second=search_sca_selector(train,validation,config);
        CHECK_TRUE(first.train_best.is_exact());
        CHECK_TRUE(first.validation_best.is_exact());
        CHECK_EQ(first.validation_best.algorithm.selector,second.validation_best.algorithm.selector);
        CHECK_TRUE(close_metric(
            first.validation_best.mean_proxy_cost,second.validation_best.mean_proxy_cost));

        const std::filesystem::path path="test_selector_roundtrip.txt";
        save_selector_genome(first.validation_best.algorithm.selector,path);
        const SelectorGenome loaded=load_selector_genome(path);
        CHECK_EQ(loaded,first.validation_best.algorithm.selector);
        std::filesystem::remove(path);
    });
}
