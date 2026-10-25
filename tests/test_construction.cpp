#include "test_util.hpp"
#include <random>
using namespace polar_scl;
int main()
{
    return run_test_main("construction",[]
    {
        CHECK_EQ(polarization_weight_order(8),(std::vector<std::size_t>
        {
            0,1,2,4,3,5,6,7
        }

        ));
        std::mt19937_64 rng(5);
        std::uniform_int_distribution<int>b(0,1);
        for(std::size_t n:{1u,7u,24u,56u})
        {
            std::vector<Bit>p(n);
            for(auto&x:p)x=static_cast<Bit>(b(rng));
            auto info=append_crc(p);
            CHECK_TRUE(crc_ok(info));
            auto bad=info;
            bad[n/2]^=1;
            CHECK_TRUE(!crc_ok(bad));
        }

        auto frozen=pw_frozen_pattern(64,40);
        std::vector<Bit>payload(32);
        for(std::size_t i=0;i<32;++i)payload[i]=static_cast<Bit>(i&1);
        auto [u,x]=encode_payload(payload,frozen);
        CHECK_EQ(u.size(),64u);
        CHECK_EQ(x.size(),64u);
        auto info=extract_information(u,frozen);
        CHECK_TRUE(std::equal(payload.begin(),payload.end(),info.begin()));
        CHECK_TRUE(crc_ok(info));
        auto f16=pw_frozen_pattern(16,12);
        std::vector<Bit>p4
        {
            1,0,1,1
        };
        auto [good,_x]=encode_payload(p4,f16);
        auto bad=good;
        auto it=std::find(f16.begin(),f16.end(),false);
        bad[static_cast<std::size_t>(std::distance(f16.begin(),it))]^=1;
        auto selected=select_crc_path({{0.0,bad},{1.0,good}},f16);
        CHECK_EQ(selected.u,good);
    }

    );
}
