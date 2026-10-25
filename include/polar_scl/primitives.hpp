#pragma once

#include "polar_scl/algorithm_ir.hpp"

#include <string>
#include <vector>

namespace polar_scl
{
    struct PrimitiveInfo
    {
        Rate1Engine engine;
        std::string name;
        std::string family;
        bool exact;
        std::string description;
    };

    inline std::vector<PrimitiveInfo> rate1_primitives()
    {
        return {
            {
                Rate1Engine::Exhaustive, "exhaustive", "oracle", true, "Enumerate every local word; small-node correctness oracle."
            }
            ,
            {
                Rate1Engine::LazyExact, "lazy_exact", "best-first", true, "Exact lazy subset-sum streams merged by a global heap."
            }
            ,
            {
                Rate1Engine::FastSSCLBounded, "fast_sscl_bounded", "Fast-SSCL-inspired", true, "Software baseline restricted to min(L-1,Ns) least-reliable positions."
            }
            ,
            {
                Rate1Engine::CertifiedFrontier, "certified_frontier", "adaptive-certified", true, "Adaptive reliability frontier with hidden-candidate certificate and exact fallback."
            }
            ,
        };
    }
}

// namespace polar_scl
