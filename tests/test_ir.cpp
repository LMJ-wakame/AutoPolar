#include "test_util.hpp"
using namespace polar_scl;
int main()
{
    return run_test_main("ir",[]
    {
        validate_policy(Rate1Policy{});
        CHECK_THROW(validate_policy(Rate1Policy{1,0,ExpansionMode::MinimumBound,64}));
        CHECK_THROW(validate_policy(Rate1Policy{1,1,ExpansionMode::MinimumBound,0}));
        validate_rate1_ir(Rate1NodeIR{});
        SelectorGenome g;
        validate_selector_genome(g);
    }

    );
}
