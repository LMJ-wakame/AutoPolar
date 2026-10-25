#include "polar_scl/algorithm_ir.hpp"

#include <stdexcept>

namespace polar_scl
{
    std::string to_string(Rate1Engine engine)
    {
        switch(engine)
        {
            case Rate1Engine::Exhaustive: return "exhaustive";
            case Rate1Engine::LazyExact: return "lazy_exact";
            case Rate1Engine::FastSSCLBounded: return "fast_sscl_bounded";
            case Rate1Engine::CertifiedFrontier: return "certified_frontier";
        }

        throw std::logic_error("unreachable Rate-1 engine");
    }

    Rate1Policy Rate1NodeIR::frontier_policy() const
    {
        Rate1Policy policy{initial_active_bits, growth_step, expansion_mode, candidate_budget};
        validate_policy(policy);
        return policy;
    }

    Rate1NodeIR SelectorGenome::ir_for_engine(Rate1Engine engine) const
    {
        return {engine, frontier_initial_active_bits, frontier_growth_step, frontier_expansion_mode, frontier_candidate_budget};
    }

    void validate_rate1_ir(const Rate1NodeIR& ir)
    {
        (void)ir.frontier_policy();
    }

    void validate_selector_genome(const SelectorGenome& genome)
    {
        validate_rate1_ir(genome.ir_for_engine(genome.low_list_engine));
        if(genome.small_node_threshold < 1) throw std::invalid_argument("small_node_threshold must be positive");
        if(genome.reliability_gap_threshold < 0.0) throw std::invalid_argument("reliability_gap_threshold must be non-negative");
        if(genome.parent_metric_gap_threshold < 0.0) throw std::invalid_argument("parent_metric_gap_threshold must be non-negative");
        validate_rate1_ir(genome.ir_for_engine(genome.small_engine));
        validate_rate1_ir(genome.ir_for_engine(genome.easy_engine));
        validate_rate1_ir(genome.ir_for_engine(genome.hard_engine));
    }

    void validate_algorithm_ir(const PolarAlgorithmIR& ir)
    {
        validate_rate1_ir(ir.fixed_rate1);
        validate_selector_genome(ir.selector);
    }
}

// namespace polar_scl
