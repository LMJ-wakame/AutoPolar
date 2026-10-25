#include "polar_scl/polar_scl.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
using namespace polar_scl;
static bool eq(const std::vector<DecoderPath>&a,const std::vector<DecoderPath>&b)
{
    if(a.size()!=b.size())return false;
    for(std::size_t i=0;i<a.size();++i)if(a[i].u!=b[i].u||std::abs(a[i].metric-b[i].metric)>1e-10*std::max(
    {
        1.0,std::abs(a[i].metric),std::abs(b[i].metric)
    }

    ))return false;
    return true;
}

int main()
{
    constexpr std::size_t n=64,payload_bits=32,crc_bits=8,k=40,L=8;
    auto frozen=pw_frozen_pattern(n,k);
    PolarAlgorithmIR algorithm;
    for(int s=0;s<3;++s)
    {
        double snr=s,sigma2=1.0/(2*std::pow(10.0,snr/10.0)),sigma=std::sqrt(sigma2);
        std::mt19937_64 rng(500+s);
        std::uniform_int_distribution<int>bd(0,1);
        std::normal_distribution<double>nd(0,1);
        std::size_t mism=0,re=0,pe=0;
        double rl=0,pl=0,rs=0,ps=0;
        std::map<Rate1Engine,std::size_t>eng;
        for(int f=0;f<50;++f)
        {
            std::vector<Bit>payload(payload_bits);
            for(auto&b:payload)b=bd(rng);
            auto [_u,x]=encode_payload(payload,frozen,crc_bits);
            std::vector<double>llrs(n);
            for(std::size_t i=0;i<n;++i)
            {
                double y=1.0-2.0*x[i]+sigma*nd(rng);
                llrs[i]=2*y/sigma2;
            }

            auto ref=scl_decode(llrs,frozen,L);
            auto port=portfolio_scl_decode(llrs,frozen,L,algorithm);
            mism+=eq(ref.paths,port.paths)?0:1;
            auto rp=select_crc_path(ref.paths,frozen,crc_bits),pp=select_crc_path(port.paths,frozen,crc_bits);
            auto ri=extract_information(rp.u,frozen),pi=extract_information(pp.u,frozen);
            re+=!std::equal(payload.begin(),payload.end(),ri.begin());
            pe+=!std::equal(payload.begin(),payload.end(),pi.begin());
            rl+=ref.stats.leaf_steps;
            pl+=port.stats.decoder.leaf_steps;
            rs+=ref.stats.sorts;
            ps+=port.stats.decoder.sorts;
            for(auto[e,c]:port.stats.engine_counts)eng[e]+=c;
        }

        std::cout<<std::fixed<<std::setprecision(1)<<"SNR="<<snr<<" dB "<<std::setprecision(3)<<"FER(ref/portfolio)="<<static_cast<double>(re)/50<<'/'<<static_cast<double>(pe)/50<<" mismatch="<<mism<<std::setprecision(1)<<" leaf="<<rl/50<<'/'<<pl/50<<" prune-events="<<rs/50<<'/'<<ps/50<<" engines=";
        for(auto[e,c]:eng)std::cout<<to_string(e)<<':'<<c<<' ';
        std::cout<<'\n';
    }
}
