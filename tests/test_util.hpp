#pragma once

#include "polar_scl/polar_scl.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define CHECK_TRUE(expr) \
    do \
    { \
        if(!(expr)) \
        { \
            throw std::runtime_error(std::string("CHECK failed: ")+#expr); \
        } \
    } while(false)

#define CHECK_EQ(a,b) \
    do \
    { \
        if(!((a)==(b))) \
        { \
            throw std::runtime_error(std::string("CHECK_EQ failed: ")+#a+" != "+#b); \
        } \
    } while(false)

#define CHECK_THROW(expr) \
    do \
    { \
        bool caught=false; \
        try \
        { \
            (void)(expr); \
        } \
        catch(...) \
        { \
            caught=true; \
        } \
        if(!caught) \
        { \
            throw std::runtime_error(std::string("expected exception: ")+#expr); \
        } \
    } while(false)

inline bool close_metric(double a,double b,double tolerance=1e-12)
{
    return std::abs(a-b)<=tolerance*std::max({1.0,std::abs(a),std::abs(b)});
}

inline void check_candidates_equal(
    const std::vector<polar_scl::Rate1Candidate>& lhs,
    const std::vector<polar_scl::Rate1Candidate>& rhs)
{
    CHECK_EQ(lhs.size(),rhs.size());
    for(std::size_t index=0;index<lhs.size();++index)
    {
        CHECK_TRUE(polar_scl::local_rank_key(lhs[index])==polar_scl::local_rank_key(rhs[index]));
    }
}

inline void check_paths_equal(
    const std::vector<polar_scl::DecoderPath>& lhs,
    const std::vector<polar_scl::DecoderPath>& rhs,
    double tolerance=1e-12)
{
    CHECK_EQ(lhs.size(),rhs.size());
    for(std::size_t index=0;index<lhs.size();++index)
    {
        CHECK_EQ(lhs[index].u,rhs[index].u);
        CHECK_TRUE(close_metric(lhs[index].metric,rhs[index].metric,tolerance));
    }
}

inline int run_test_main(const char* name,const std::function<void()>& function)
{
    try
    {
        function();
        std::cout<<"[PASS] "<<name<<'\n';
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"[FAIL] "<<name<<": "<<error.what()<<'\n';
        return 1;
    }
}
