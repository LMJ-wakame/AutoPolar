#include "polar_scl/stress.hpp"

#include "polar_scl/construction.hpp"
#include "polar_scl/portfolio.hpp"
#include "polar_scl/rate1.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace polar_scl
{
    namespace
    {
        template<typename T>
        void write_binary(std::ostream& output,const T& value)
        {
            static_assert(std::is_trivially_copyable_v<T>);
            output.write(reinterpret_cast<const char*>(&value),static_cast<std::streamsize>(sizeof(T)));
            if(!output)
            {
                throw std::runtime_error("failed to write stress dataset");
            }
        }

        template<typename T>
        T read_binary(std::istream& input)
        {
            static_assert(std::is_trivially_copyable_v<T>);
            T value{};
            input.read(reinterpret_cast<char*>(&value),static_cast<std::streamsize>(sizeof(T)));
            if(!input)
            {
                throw std::runtime_error("unexpected end of stress dataset");
            }
            return value;
        }

        double mean(const std::vector<double>& values)
        {
            if(values.empty())
            {
                return 0.0;
            }
            return std::accumulate(values.begin(),values.end(),0.0)/static_cast<double>(values.size());
        }

        Rate1NodeIR fixed_rate1_ir(Rate1Engine engine)
        {
            Rate1NodeIR ir;
            ir.engine=engine;
            ir.initial_active_bits=0;
            ir.growth_step=1;
            ir.expansion_mode=ExpansionMode::MinimumBound;
            ir.candidate_budget=256;
            return ir;
        }

        PolarAlgorithmIR algorithm_for_variant(
            DecoderVariant variant,
            const PolarAlgorithmIR& searched_algorithm)
        {
            if(variant==DecoderVariant::SearchedPortfolio)
            {
                return searched_algorithm;
            }
            if(variant==DecoderVariant::PortfolioRule)
            {
                return PolarAlgorithmIR{};
            }

            PolarAlgorithmIR algorithm;
            algorithm.selector_mode=SelectorMode::Fixed;
            switch(variant)
            {
                case DecoderVariant::SSCLRate01Exact:
                algorithm.fixed_rate1=fixed_rate1_ir(Rate1Engine::LazyExact);
                return algorithm;
                case DecoderVariant::FastSSCLRate1Bounded:
                algorithm.fixed_rate1=fixed_rate1_ir(Rate1Engine::FastSSCLBounded);
                return algorithm;
                case DecoderVariant::CertifiedFrontier:
                algorithm.fixed_rate1=fixed_rate1_ir(Rate1Engine::CertifiedFrontier);
                return algorithm;
                default:
                throw std::invalid_argument("stress direct mode requires a fast portfolio decoder");
            }
        }

        std::vector<bool> make_frozen(
            StressMode mode,
            std::size_t length,
            std::mt19937_64& rng)
        {
            if(mode==StressMode::RepetitionWholeNode)
            {
                std::vector<bool> frozen(length,true);
                frozen.back()=false;
                return frozen;
            }
            if(mode==StressMode::SpcWholeNode)
            {
                std::vector<bool> frozen(length,false);
                frozen.front()=true;
                return frozen;
            }
            if(mode==StressMode::ParentCompetition)
            {
                std::vector<bool> frozen(length,false);
                const std::size_t prefix=std::max<std::size_t>(4,length/4);
                for(std::size_t index=0;index<prefix;++index)
                {
                    frozen[index]=(index%3)==0;
                }
                return frozen;
            }
            if(mode==StressMode::RandomFrozen||mode==StressMode::DeepRandomFrozen)
            {
                const std::size_t low=mode==StressMode::DeepRandomFrozen?length/3:1;
                const std::size_t high=mode==StressMode::DeepRandomFrozen?2*length/3:length-1;
                std::uniform_int_distribution<std::size_t> count_distribution(
                    std::max<std::size_t>(1,low),
                    std::max<std::size_t>(1,high));
                const std::size_t information_bits=count_distribution(rng);
                std::vector<std::size_t> indices(length);
                std::iota(indices.begin(),indices.end(),std::size_t{0});
                std::shuffle(indices.begin(),indices.end(),rng);
                std::vector<bool> frozen(length,true);
                for(std::size_t index=0;index<information_bits;++index)
                {
                    frozen[indices[index]]=false;
                }
                return frozen;
            }

            const std::array<double,3> rates{0.25,0.5,0.75};
            std::uniform_int_distribution<std::size_t> distribution(0,rates.size()-1);
            const auto information_bits=static_cast<std::size_t>(
                std::llround(rates[distribution(rng)]*static_cast<double>(length)));
            return pw_frozen_pattern(length,std::max<std::size_t>(1,information_bits));
        }

        double quantize(double value,int bits)
        {
            const int limit=(1<<(bits-1))-1;
            const double scaled=std::round(value*2.0);
            return std::clamp(scaled,static_cast<double>(-limit),static_cast<double>(limit));
        }

        double random_sign(std::mt19937_64& rng)
        {
            return (rng()&1ULL)==0ULL?-1.0:1.0;
        }

        std::vector<double> make_llrs(
            StressMode mode,
            std::size_t length,
            std::mt19937_64& rng)
        {
            std::normal_distribution<double> normal(0.0,3.0);
            std::uniform_real_distribution<double> tiny(-0.02,0.02);
            std::uniform_int_distribution<int> tie_value(-2,2);
            std::vector<double> llrs(length);

            for(std::size_t index=0;index<length;++index)
            {
                const double base=normal(rng);
                switch(mode)
                {
                    case StressMode::Continuous:
                    case StressMode::RepetitionWholeNode:
                    case StressMode::SpcWholeNode:
                    case StressMode::RandomFrozen:
                    llrs[index]=base;
                    break;
                    case StressMode::Quantized8:
                    llrs[index]=quantize(base,8);
                    break;
                    case StressMode::Quantized6:
                    llrs[index]=quantize(base,6);
                    break;
                    case StressMode::Quantized4:
                    llrs[index]=quantize(base,4);
                    break;
                    case StressMode::Quantized3:
                    llrs[index]=quantize(base,3);
                    break;
                    case StressMode::NearZero:
                    llrs[index]=tiny(rng);
                    break;
                    case StressMode::TieHeavy:
                    llrs[index]=static_cast<double>(tie_value(rng));
                    break;
                    case StressMode::MetricBoundary:
                    {
                        static constexpr double levels[]{0.0,0.25,-0.25,0.5,-0.5,1.0,-1.0};
                        llrs[index]=levels[index%(sizeof(levels)/sizeof(levels[0]))];
                        if(index%11==0)
                        {
                            llrs[index]+=static_cast<double>(static_cast<int>(index%3)-1)*1e-12;
                        }
                        break;
                    }
                    case StressMode::SubsetSumCollision:
                    llrs[index]=random_sign(rng)*static_cast<double>(index+1);
                    break;
                    case StressMode::GuardBoundaryOutside:
                    llrs[index]=random_sign(rng)*(1.0+static_cast<double>(index)*1.35e-10);
                    break;
                    case StressMode::ParentCompetition:
                    llrs[index]=random_sign(rng)*(0.08+0.017*static_cast<double>(index+1));
                    break;
                    case StressMode::DynamicRange:
                    {
                        const double exponent=-8.0+16.0*static_cast<double>(index)/
                            static_cast<double>(std::max<std::size_t>(1,length-1));
                        llrs[index]=random_sign(rng)*std::pow(10.0,exponent);
                        break;
                    }
                    case StressMode::SignedZero:
                    llrs[index]=(index&1U)==0U?0.0:-0.0;
                    break;
                    case StressMode::CancellationChain:
                    {
                        const double magnitude=1.0+0.003*static_cast<double>(index+1);
                        llrs[index]=(index&1U)==0U?magnitude:-(magnitude+2.2e-10);
                        break;
                    }
                    case StressMode::DeepRandomFrozen:
                    llrs[index]=random_sign(rng)*(0.02+0.011*static_cast<double>(index+1));
                    break;
                    case StressMode::PermutedReliability:
                    {
                        const double magnitude=0.125+0.03125*static_cast<double>(2*index+1);
                        llrs[index]=random_sign(rng)*magnitude;
                        break;
                    }
                }
            }

            if(mode!=StressMode::MetricBoundary&&mode!=StressMode::SignedZero)
            {
                std::shuffle(llrs.begin(),llrs.end(),rng);
            }
            return llrs;
        }

        PortfolioDecodeResult direct_decode(
            DecoderVariant variant,
            const StressCase& test_case,
            const PolarAlgorithmIR& searched_algorithm,
            bool raw)
        {
            const PolarAlgorithmIR algorithm=algorithm_for_variant(variant,searched_algorithm);
            PortfolioDecodeOptions options;
            if(raw)
            {
                options.enable_rate1_tie_guard=false;
                options.enable_spc_cutoff_guard=false;
                options.enable_final_list_tie_guard=false;
            }
            return portfolio_scl_decode(
                test_case.llrs,
                test_case.frozen,
                test_case.list_size,
                algorithm,
                options);
        }
    }

    std::string to_string(StressMode mode)
    {
        switch(mode)
        {
            case StressMode::Continuous:
            return "continuous";
            case StressMode::Quantized8:
            return "quantized_8bit";
            case StressMode::Quantized6:
            return "quantized_6bit";
            case StressMode::Quantized4:
            return "quantized_4bit";
            case StressMode::Quantized3:
            return "quantized_3bit";
            case StressMode::NearZero:
            return "near_zero";
            case StressMode::TieHeavy:
            return "tie_heavy";
            case StressMode::MetricBoundary:
            return "metric_boundary";
            case StressMode::RepetitionWholeNode:
            return "repetition_whole_node";
            case StressMode::SpcWholeNode:
            return "spc_whole_node";
            case StressMode::RandomFrozen:
            return "random_frozen";
            case StressMode::SubsetSumCollision:
            return "subset_sum_collision";
            case StressMode::GuardBoundaryOutside:
            return "guard_boundary_outside";
            case StressMode::ParentCompetition:
            return "parent_competition";
            case StressMode::DynamicRange:
            return "dynamic_range";
            case StressMode::SignedZero:
            return "signed_zero";
            case StressMode::CancellationChain:
            return "cancellation_chain";
            case StressMode::DeepRandomFrozen:
            return "deep_random_frozen";
            case StressMode::PermutedReliability:
            return "permuted_reliability";
        }
        throw std::logic_error("unreachable stress mode");
    }

    std::string to_string(StressLayer layer)
    {
        switch(layer)
        {
            case StressLayer::GuardedProduction:
            return "guarded";
            case StressLayer::DirectInnerGuard:
            return "direct";
            case StressLayer::RawFastPath:
            return "raw";
        }
        throw std::logic_error("unreachable stress layer");
    }

    StressDataset generate_stress_dataset(const StressDatasetConfig& config)
    {
        if(config.cases_per_mode==0||config.lengths.empty()||config.list_sizes.empty())
        {
            throw std::invalid_argument("stress configuration must be non-empty");
        }

        const std::vector<StressMode> modes{
            StressMode::Continuous,
            StressMode::Quantized8,
            StressMode::Quantized6,
            StressMode::Quantized4,
            StressMode::Quantized3,
            StressMode::NearZero,
            StressMode::TieHeavy,
            StressMode::MetricBoundary,
            StressMode::RepetitionWholeNode,
            StressMode::SpcWholeNode,
            StressMode::RandomFrozen,
            StressMode::SubsetSumCollision,
            StressMode::GuardBoundaryOutside,
            StressMode::ParentCompetition,
            StressMode::DynamicRange,
            StressMode::SignedZero,
            StressMode::CancellationChain,
            StressMode::DeepRandomFrozen,
            StressMode::PermutedReliability};

        std::mt19937_64 rng(config.seed);
        std::uniform_int_distribution<std::size_t> length_distribution(0,config.lengths.size()-1);
        std::uniform_int_distribution<std::size_t> list_distribution(0,config.list_sizes.size()-1);

        StressDataset dataset;
        dataset.seed=config.seed;
        dataset.cases.reserve(config.cases_per_mode*modes.size());
        for(StressMode mode:modes)
        {
            for(std::size_t case_index=0;case_index<config.cases_per_mode;++case_index)
            {
                std::size_t length=config.lengths[length_distribution(rng)];
                if((mode==StressMode::RepetitionWholeNode||mode==StressMode::SpcWholeNode)&&length>64)
                {
                    length=64;
                }
                if(length<2||(length&(length-1))!=0||length>128)
                {
                    throw std::invalid_argument("stress lengths must be powers of two in [2,128]");
                }
                StressCase test_case;
                test_case.mode=mode;
                test_case.list_size=std::min(config.list_sizes[list_distribution(rng)],length);
                test_case.frozen=make_frozen(mode,length,rng);
                test_case.llrs=make_llrs(mode,length,rng);
                dataset.cases.push_back(std::move(test_case));
            }
        }
        return dataset;
    }

    void save_stress_dataset_binary(const StressDataset& dataset,const std::filesystem::path& path)
    {
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }
        std::ofstream output(path,std::ios::binary);
        if(!output)
        {
            throw std::runtime_error("failed to create stress dataset: "+path.string());
        }

        constexpr char magic[8]={'S','C','A','S','T','R','0','2'};
        output.write(magic,8);
        write_binary(output,dataset.seed);
        write_binary(output,static_cast<std::uint64_t>(dataset.cases.size()));
        for(const StressCase& test_case:dataset.cases)
        {
            if(test_case.frozen.size()!=test_case.llrs.size()||test_case.llrs.empty())
            {
                throw std::invalid_argument("invalid stress case shape");
            }
            write_binary(output,static_cast<std::uint32_t>(test_case.mode));
            write_binary(output,static_cast<std::uint32_t>(test_case.list_size));
            write_binary(output,static_cast<std::uint32_t>(test_case.llrs.size()));
            for(bool frozen:test_case.frozen)
            {
                write_binary(output,static_cast<std::uint8_t>(frozen?1:0));
            }
            output.write(
                reinterpret_cast<const char*>(test_case.llrs.data()),
                static_cast<std::streamsize>(test_case.llrs.size()*sizeof(double)));
            if(!output)
            {
                throw std::runtime_error("failed to write stress LLRs");
            }
        }
    }

    StressDataset load_stress_dataset_binary(const std::filesystem::path& path)
    {
        std::ifstream input(path,std::ios::binary);
        if(!input)
        {
            throw std::runtime_error("failed to open stress dataset: "+path.string());
        }
        char magic[8]{};
        input.read(magic,8);
        constexpr char expected[8]={'S','C','A','S','T','R','0','2'};
        if(!input||std::memcmp(magic,expected,8)!=0)
        {
            throw std::runtime_error("invalid stress dataset magic/version");
        }

        StressDataset dataset;
        dataset.seed=read_binary<std::uint64_t>(input);
        const auto count=static_cast<std::size_t>(read_binary<std::uint64_t>(input));
        dataset.cases.reserve(count);
        for(std::size_t case_index=0;case_index<count;++case_index)
        {
            StressCase test_case;
            test_case.mode=static_cast<StressMode>(read_binary<std::uint32_t>(input));
            test_case.list_size=static_cast<std::size_t>(read_binary<std::uint32_t>(input));
            const auto length=static_cast<std::size_t>(read_binary<std::uint32_t>(input));
            test_case.frozen.resize(length);
            for(std::size_t index=0;index<length;++index)
            {
                test_case.frozen[index]=read_binary<std::uint8_t>(input)!=0;
            }
            test_case.llrs.resize(length);
            input.read(
                reinterpret_cast<char*>(test_case.llrs.data()),
                static_cast<std::streamsize>(length*sizeof(double)));
            if(!input)
            {
                throw std::runtime_error("unexpected end of stress LLRs");
            }
            dataset.cases.push_back(std::move(test_case));
        }
        return dataset;
    }

    double StressValidationRow::guarded_mismatch_rate() const
    {
        return cases==0?0.0:static_cast<double>(guarded_mismatches)/static_cast<double>(cases);
    }

    double StressValidationRow::direct_mismatch_rate() const
    {
        return cases==0?0.0:static_cast<double>(direct_mismatches)/static_cast<double>(cases);
    }

    double StressValidationRow::raw_mismatch_rate() const
    {
        return cases==0?0.0:static_cast<double>(raw_mismatches)/static_cast<double>(cases);
    }

    std::size_t StressValidationRow::production_fast_path_cases() const
    {
        return cases>=outer_guard_cases?cases-outer_guard_cases:0;
    }

    StressValidationReport run_stress_validation(
        const StressDataset& dataset,
        const PolarAlgorithmIR& searched_algorithm,
        bool include_searched,
        std::size_t max_saved_mismatches)
    {
        std::vector<DecoderVariant> variants{
            DecoderVariant::SSCLRate01Exact,
            DecoderVariant::FastSSCLRate1Bounded,
            DecoderVariant::CertifiedFrontier,
            DecoderVariant::PortfolioRule};
        if(include_searched)
        {
            variants.push_back(DecoderVariant::SearchedPortfolio);
        }

        struct Accumulator
        {
            StressValidationRow row;
            std::vector<double> guarded_times;
            std::vector<double> direct_times;
            std::vector<double> raw_times;
        };
        std::map<std::pair<int,int>,Accumulator> accumulators;
        StressValidationReport report;

        for(const StressCase& test_case:dataset.cases)
        {
            const auto reference=scl_decode(test_case.llrs,test_case.frozen,test_case.list_size);
            const ParentPath probe{0,0.0,test_case.llrs};
            const bool outer_guard=rate1_metric_tie_risk({probe});

            for(DecoderVariant variant:variants)
            {
                const auto key=std::make_pair(static_cast<int>(variant),static_cast<int>(test_case.mode));
                auto& accumulator=accumulators[key];
                accumulator.row.variant=variant;
                accumulator.row.mode=test_case.mode;
                ++accumulator.row.cases;
                accumulator.row.outer_guard_cases+=outer_guard?1U:0U;

                const PolarAlgorithmIR algorithm=variant==DecoderVariant::SearchedPortfolio
                    ?searched_algorithm:PolarAlgorithmIR{};

                const auto guarded_begin=std::chrono::steady_clock::now();
                const auto guarded=decode_variant(
                    variant,
                    test_case.llrs,
                    test_case.frozen,
                    test_case.list_size,
                    algorithm);
                const auto guarded_end=std::chrono::steady_clock::now();
                accumulator.guarded_times.push_back(std::chrono::duration<double,std::micro>(
                    guarded_end-guarded_begin).count());
                const bool guarded_equal=paths_equivalent(reference.paths,guarded.paths,1e-10);
                accumulator.row.guarded_mismatches+=guarded_equal?0U:1U;
                if(!guarded_equal&&report.mismatch_examples.size()<max_saved_mismatches)
                {
                    report.mismatch_examples.push_back({
                        variant,
                        test_case.mode,
                        StressLayer::GuardedProduction,
                        test_case.list_size,
                        test_case.frozen,
                        test_case.llrs});
                }
                accumulator.row.guarded_tie_fallbacks+=guarded.tie_fallbacks;
                accumulator.row.guarded_rate1_fallbacks+=guarded.stats.rate1_fallbacks;

                const auto direct_begin=std::chrono::steady_clock::now();
                const auto direct=direct_decode(variant,test_case,searched_algorithm,false);
                const auto direct_end=std::chrono::steady_clock::now();
                accumulator.direct_times.push_back(std::chrono::duration<double,std::micro>(
                    direct_end-direct_begin).count());
                const bool direct_equal=paths_equivalent(reference.paths,direct.paths,1e-10);
                accumulator.row.direct_mismatches+=direct_equal?0U:1U;
                if(!direct_equal&&report.mismatch_examples.size()<max_saved_mismatches)
                {
                    report.mismatch_examples.push_back({
                        variant,
                        test_case.mode,
                        StressLayer::DirectInnerGuard,
                        test_case.list_size,
                        test_case.frozen,
                        test_case.llrs});
                }
                accumulator.row.direct_tie_fallbacks+=direct.stats.tie_fallbacks;
                accumulator.row.direct_rate1_fallbacks+=direct.stats.decoder.rate1_fallbacks;

                const auto raw_begin=std::chrono::steady_clock::now();
                const auto raw=direct_decode(variant,test_case,searched_algorithm,true);
                const auto raw_end=std::chrono::steady_clock::now();
                accumulator.raw_times.push_back(std::chrono::duration<double,std::micro>(
                    raw_end-raw_begin).count());
                const bool raw_equal=paths_equivalent(reference.paths,raw.paths,1e-10);
                accumulator.row.raw_mismatches+=raw_equal?0U:1U;
                if(!raw_equal&&report.mismatch_examples.size()<max_saved_mismatches)
                {
                    report.mismatch_examples.push_back({
                        variant,
                        test_case.mode,
                        StressLayer::RawFastPath,
                        test_case.list_size,
                        test_case.frozen,
                        test_case.llrs});
                }
            }
        }

        report.rows.reserve(accumulators.size());
        for(auto& [key,accumulator]:accumulators)
        {
            (void)key;
            accumulator.row.guarded_mean_us=mean(accumulator.guarded_times);
            accumulator.row.direct_mean_us=mean(accumulator.direct_times);
            accumulator.row.raw_mean_us=mean(accumulator.raw_times);
            report.rows.push_back(accumulator.row);
        }
        return report;
    }

    void write_stress_validation_csv(
        const std::vector<StressValidationRow>& rows,
        const std::filesystem::path& path)
    {
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }
        std::ofstream output(path);
        if(!output)
        {
            throw std::runtime_error("failed to create stress validation CSV");
        }
        output<<"decoder,mode,cases,outer_guard_cases,production_fast_path_cases,"
            <<"guarded_mismatches,guarded_mismatch_rate,guarded_tie_fallbacks,guarded_rate1_fallbacks,guarded_mean_us,"
            <<"direct_mismatches,direct_mismatch_rate,direct_tie_fallbacks,direct_rate1_fallbacks,direct_mean_us,"
            <<"raw_mismatches,raw_mismatch_rate,raw_mean_us\n";
        output<<std::setprecision(12);
        for(const StressValidationRow& row:rows)
        {
            output<<to_string(row.variant)<<','
                <<to_string(row.mode)<<','
                <<row.cases<<','
                <<row.outer_guard_cases<<','
                <<row.production_fast_path_cases()<<','
                <<row.guarded_mismatches<<','
                <<row.guarded_mismatch_rate()<<','
                <<row.guarded_tie_fallbacks<<','
                <<row.guarded_rate1_fallbacks<<','
                <<row.guarded_mean_us<<','
                <<row.direct_mismatches<<','
                <<row.direct_mismatch_rate()<<','
                <<row.direct_tie_fallbacks<<','
                <<row.direct_rate1_fallbacks<<','
                <<row.direct_mean_us<<','
                <<row.raw_mismatches<<','
                <<row.raw_mismatch_rate()<<','
                <<row.raw_mean_us<<'\n';
        }
    }

    void write_stress_mismatch_examples_csv(
        const std::vector<StressMismatchExample>& examples,
        const std::filesystem::path& path)
    {
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }
        std::ofstream output(path);
        if(!output)
        {
            throw std::runtime_error("failed to create mismatch example CSV");
        }
        output<<"decoder,mode,layer,list_size,length,frozen,llrs\n";
        output<<std::setprecision(17);
        for(const auto& example:examples)
        {
            std::string frozen;
            frozen.reserve(example.frozen.size());
            for(bool value:example.frozen)
            {
                frozen.push_back(value?'1':'0');
            }
            std::ostringstream llrs;
            for(std::size_t index=0;index<example.llrs.size();++index)
            {
                if(index>0)
                {
                    llrs<<';';
                }
                llrs<<std::setprecision(17)<<example.llrs[index];
            }
            output<<to_string(example.variant)<<','
                <<to_string(example.mode)<<','
                <<to_string(example.layer)<<','
                <<example.list_size<<','
                <<example.llrs.size()<<','
                <<frozen<<','
                <<'"'<<llrs.str()<<'"'<<'\n';
        }
    }
}
