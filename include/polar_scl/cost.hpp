#pragma once

#include "polar_scl/rate1.hpp"

#include <cstddef>

namespace polar_scl
{
    struct CostWeights
    {
        double candidate{1.0};
        double heap_operation{0.20};
        double frontier_expansion{0.50};
        double certificate_check{0.10};
    };

    struct ComplexityVector
    {
        std::size_t candidates{};
        std::size_t heap_pushes{};
        std::size_t heap_pops{};
        std::size_t frontier_expansions{};
        std::size_t certificate_checks{};
        bool fallback_used{false};

        std::size_t heap_operations() const
        {
            return heap_pushes + heap_pops;
        }

        double proxy_cost(const CostWeights& weights = {}) const;
    };

    ComplexityVector complexity_from_stats(const SearchStats& stats);
}

// namespace polar_scl
