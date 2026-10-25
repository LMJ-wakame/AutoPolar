#include "polar_scl/polar_scl.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
using namespace polar_scl;
struct Frame
{
    std::vector<double>llrs;
};
static double percentile(std::vector<double>v,double q)
{
    std::sort(v.begin(),v.end());
    auto i=static_cast<std::size_t>(q*static_cast<double>(v.size()-1));
    return v[i];
}

template<class F>static void bench(const char*name,const std::vector<Frame>&frames,F fn)
{
    std::vector<double>us;
    us.reserve(frames.size());
    volatile double checksum=0;
    for(const auto&f:frames)
    {
        auto t0=std::chrono::steady_clock::now();
        auto r=fn(f.llrs);
        auto t1=std::chrono::steady_clock::now();
        checksum+=r.paths.front().metric;
        us.push_back(std::chrono::duration<double,std::micro>(t1-t0).count());
    }

    double mean=std::accumulate(us.begin(),us.end(),0.0)/us.size();
    std::cout<<std::left<<std::setw(22)<<name<<std::right<<" mean_us/frame="<<std::fixed<<std::setprecision(2)<<mean<<" p50="<<percentile(us,.50)<<" p95="<<percentile(us,.95)<<" p99="<<percentile(us,.99)<<" checksum="<<checksum<<'\n';
}

int main()
{
    constexpr std::size_t n=64,k=40,L=8,frames_n=120;
    auto frozen=pw_frozen_pattern(n,k);
    double snr=1.0,sigma2=1.0/(2*std::pow(10.0,snr/10.0)),sigma=std::sqrt(sigma2);
    std::mt19937_64 rng(9001);
    std::uniform_int_distribution<int>bd(0,1);
    std::normal_distribution<double>nd(0,1);
    std::vector<Frame>frames;
    for(std::size_t f=0;f<frames_n;++f)
    {
        std::vector<Bit>payload(32);
        for(auto&b:payload)b=bd(rng);
        auto[_u,x]=encode_payload(payload,frozen,8);
        Frame fr;
        fr.llrs.resize(n);
        for(std::size_t i=0;i<n;++i)
        {
            double y=1.0-2.0*x[i]+sigma*nd(rng);
            fr.llrs[i]=2*y/sigma2;
        }

        frames.push_back(std::move(fr));
    }

    Rate1Policy policy{0,1,ExpansionMode::MinimumBound,256};
    PolarAlgorithmIR alg;
    std::cout<<"Release benchmark, N=64 K=40 L=8 SNR=1 dB, "<<frames_n<<" pre-generated frames\n";
    bench("bitwise_scl",frames,[&](const auto&x)
    {
        return scl_decode(x,frozen,L);
    }

    );
    bench("certified_rate1",frames,[&](const auto&x)
    {
        return rate1_scl_decode(x,frozen,L,policy);
    }

    );
    bench("portfolio_scl",frames,[&](const auto&x)
    {
        auto p=portfolio_scl_decode(x,frozen,L,alg);
        return DecodeResult
        {
            p.paths,p.stats.decoder
        };
    }

    );
}
