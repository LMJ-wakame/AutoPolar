#include "test_util.hpp"

using namespace polar_scl;

int main()
{
    return run_test_main("searched_portfolio",[]
    {
        SelectorGenome genome;
        genome.small_node_threshold=4;
        genome.small_engine=Rate1Engine::LazyExact;
        genome.easy_engine=Rate1Engine::CertifiedFrontier;
        genome.hard_engine=Rate1Engine::FastSSCLBounded;
        genome.reliability_gap_threshold=0.5;
        genome.parent_metric_gap_threshold=0.25;
        genome.frontier_initial_active_bits=0;
        genome.frontier_growth_step=1;
        genome.frontier_expansion_mode=ExpansionMode::MinimumBound;
        genome.frontier_candidate_budget=256;

        PolarAlgorithmIR algorithm;
        algorithm.selector_mode=SelectorMode::Rule;
        algorithm.selector=genome;

        AwgnDatasetConfig config;
        config.length=64;
        config.payload_bits=24;
        config.crc_bits=8;
        config.ebn0_db=1.5;
        config.frames=8;
        config.seed=404;
        const Dataset dataset=generate_awgn_dataset(config);

        ExperimentRunConfig run;
        run.list_size=4;
        run.compare_with_scl=true;
        run.warmup_frames=1;
        const ExperimentMetrics metrics=run_decoder_experiment(
            dataset,DecoderVariant::SearchedPortfolio,run,algorithm);
        CHECK_EQ(metrics.list_mismatches,0U);
        CHECK_EQ(metrics.compared_frames,dataset.frames.size());
        CHECK_TRUE(metrics.mean_us>0.0);
        CHECK_TRUE(
            metrics.mean_lazy_exact_nodes+metrics.mean_fast_sscl_nodes+
            metrics.mean_certified_frontier_nodes>0.0);
    });
}
