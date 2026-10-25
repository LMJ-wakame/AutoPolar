#include "test_util.hpp"

#include <cmath>
#include <random>

using namespace polar_scl;

int main()
{
    return run_test_main("baselines",[]
    {
        std::mt19937_64 rng(20260816);
        std::normal_distribution<double> continuous(0.0,3.0);
        std::uniform_int_distribution<int> quantized(-10,10);
        const std::vector<DecoderVariant> variants{
            DecoderVariant::SSCLRate01Exact,
            DecoderVariant::FastSSCLRate1Bounded,
            DecoderVariant::CertifiedFrontier,
            DecoderVariant::PortfolioRule
        };

        for(bool use_quantized:{false,true})
        {
            for(std::size_t length:{16U,32U,64U})
            {
                for(double rate:{0.5,0.75})
                {
                    const std::size_t information_bits=static_cast<std::size_t>(
                    std::llround(rate*static_cast<double>(length)));
                    const auto frozen=pw_frozen_pattern(length,information_bits);
                    for(std::size_t list_size:{2U,4U,8U})
                    {
                        for(int trial=0;trial<18;++trial)
                        {
                            std::vector<double> llrs(length);
                            for(double& llr:llrs)
                            {
                                llr=use_quantized?static_cast<double>(quantized(rng)):continuous(rng);
                            }

                            const auto reference=scl_decode(llrs,frozen,list_size);
                            for(DecoderVariant variant:variants)
                            {
                                const auto actual=decode_variant(variant,llrs,frozen,list_size);
                                check_paths_equal(actual.paths,reference.paths,use_quantized?0.0:1e-10);
                            }
                        }
                    }
                }
            }
        }
    }
    );
}
