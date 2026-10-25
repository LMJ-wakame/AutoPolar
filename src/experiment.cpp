#include "polar_scl/experiment.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace polar_scl
{
    namespace
    {
        std::vector<Bit> random_payload(std::size_t size,std::mt19937_64& rng)
        {
            std::uniform_int_distribution<int> distribution(0,1);
            std::vector<Bit> payload(size);
            for(Bit& bit:payload)
            {
                bit=static_cast<Bit>(distribution(rng));
            }

            return payload;
        }

        std::string bits_to_string(const std::vector<Bit>& bits)
        {
            std::string result;
            result.reserve(bits.size());
            for(Bit bit:bits)
            {
                result.push_back(bit==0?'0':'1');
            }

            return result;
        }

        std::vector<Bit> string_to_bits(const std::string& text)
        {
            std::vector<Bit> result;
            result.reserve(text.size());
            for(char value:text)
            {
                if(value=='0')
                {
                    result.push_back(Bit{0});
                }

                else if(value=='1')
                {
                    result.push_back(Bit{1});
                }

                else
                {
                    throw std::invalid_argument("payload contains a non-binary character");
                }
            }

            return result;
        }

        std::vector<std::string> split_csv_line(const std::string& line)
        {
            std::vector<std::string> fields;
            std::stringstream stream(line);
            std::string field;
            while(std::getline(stream,field,','))
            {
                fields.push_back(field);
            }

            return fields;
        }

        double percentile(std::vector<double> values,double probability)
        {
            if(values.empty())
            {
                return 0.0;
            }

            std::sort(values.begin(),values.end());
            const std::size_t index=static_cast<std::size_t>(
            std::floor(probability*static_cast<double>(values.size()-1)));
            return values[index];
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

        PolarAlgorithmIR fixed_algorithm(Rate1Engine engine)
        {
            PolarAlgorithmIR algorithm;
            algorithm.selector_mode=SelectorMode::Fixed;
            algorithm.fixed_rate1=fixed_rate1_ir(engine);
            return algorithm;
        }

        std::size_t count_payload_errors(const std::vector<Bit>& expected,const std::vector<Bit>& actual)
        {
            if(actual.size()<expected.size())
            {
                return expected.size();
            }

            std::size_t errors=0;
            for(std::size_t index=0;index<expected.size();++index)
            {
                errors+=expected[index]!=actual[index]?1U:0U;
            }

            return errors;
        }

        template<typename T>
        void write_binary(std::ostream& output,const T& value)
        {
            static_assert(std::is_trivially_copyable_v<T>);
            output.write(reinterpret_cast<const char*>(&value),static_cast<std::streamsize>(sizeof(T)));
            if(!output)
            {
                throw std::runtime_error("failed to write binary dataset");
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
                throw std::runtime_error("unexpected end of binary dataset");
            }
            return value;
        }
    }
    // namespace

    std::string to_string(DecoderVariant variant)
    {
        switch(variant)
        {
            case DecoderVariant::SC:
            return "SC";
            case DecoderVariant::SCL:
            return "SCL";
            case DecoderVariant::SSCLRate01Exact:
            return "SSCL-Rate01-Exact";
            case DecoderVariant::FastSSCLRate1Bounded:
            return "FastSSCL-Rate1-Bounded";
            case DecoderVariant::CertifiedFrontier:
            return "Certified-Frontier";
            case DecoderVariant::PortfolioRule:
            return "Portfolio-Rule";
            case DecoderVariant::SearchedPortfolio:
            return "SCA-Polar";
        }

        throw std::logic_error("unreachable decoder variant");
    }

    std::vector<DecoderVariant> paper_decoder_variants(bool include_searched)
    {
        std::vector<DecoderVariant> variants{
            DecoderVariant::SC,
            DecoderVariant::SCL,
            DecoderVariant::SSCLRate01Exact,
            DecoderVariant::FastSSCLRate1Bounded,
            DecoderVariant::CertifiedFrontier,
            DecoderVariant::PortfolioRule
        };
        if(include_searched)
        {
            variants.push_back(DecoderVariant::SearchedPortfolio);
        }

        return variants;
    }

    std::size_t AwgnDatasetConfig::information_bits() const
    {
        return payload_bits+crc_bits;
    }

    double AwgnDatasetConfig::code_rate() const
    {
        return length==0?0.0:static_cast<double>(information_bits())/static_cast<double>(length);
    }

    Dataset generate_awgn_dataset(const AwgnDatasetConfig& config)
    {
        if(config.length==0||(config.length&(config.length-1))!=0)
        {
            throw std::invalid_argument("dataset length must be a positive power of two");
        }

        if(config.information_bits()==0||config.information_bits()>config.length)
        {
            throw std::invalid_argument("payload_bits+crc_bits must lie in [1,length]");
        }

        if(config.frames==0)
        {
            throw std::invalid_argument("dataset must contain at least one frame");
        }

        const auto frozen=pw_frozen_pattern(config.length,config.information_bits());
        const double ebn0_linear=std::pow(10.0,config.ebn0_db/10.0);
        const double sigma2=1.0/(2.0*config.code_rate()*ebn0_linear);
        const double sigma=std::sqrt(sigma2);
        std::mt19937_64 rng(config.seed);
        std::normal_distribution<double> noise(0.0,1.0);

        Dataset dataset;
        dataset.config=config;
        dataset.frames.reserve(config.frames);
        for(std::size_t frame_index=0;frame_index<config.frames;++frame_index)
        {
            ExperimentFrame frame;
            frame.payload=random_payload(config.payload_bits,rng);
            const auto encoded=encode_payload(frame.payload,frozen,config.crc_bits);
            const auto& codeword=encoded.second;
            frame.llrs.resize(config.length);
            for(std::size_t index=0;index<config.length;++index)
            {
                const double bpsk=1.0-2.0*static_cast<double>(codeword[index]);
                const double received=bpsk+sigma*noise(rng);
                frame.llrs[index]=2.0*received/sigma2;
            }

            dataset.frames.push_back(std::move(frame));
        }

        return dataset;
    }

    void save_dataset_csv(const Dataset& dataset,const std::filesystem::path& path)
    {
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }

        std::ofstream output(path);
        if(!output)
        {
            throw std::runtime_error("failed to create dataset file: "+path.string());
        }

        output<<"frame,length,payload_bits,crc_bits,ebn0_db,seed,payload";
        for(std::size_t index=0;index<dataset.config.length;++index)
        {
            output<<",llr_"<<index;
        }

        output<<'\n';
        output<<std::setprecision(17);
        for(std::size_t frame_index=0;frame_index<dataset.frames.size();++frame_index)
        {
            const auto& frame=dataset.frames[frame_index];
            output<<frame_index<<','
            <<dataset.config.length<<','
            <<dataset.config.payload_bits<<','
            <<dataset.config.crc_bits<<','
            <<dataset.config.ebn0_db<<','
            <<dataset.config.seed<<','
            <<bits_to_string(frame.payload);
            for(double llr:frame.llrs)
            {
                output<<','<<llr;
            }

            output<<'\n';
        }
    }

    Dataset load_dataset_csv(const std::filesystem::path& path)
    {
        std::ifstream input(path);
        if(!input)
        {
            throw std::runtime_error("failed to open dataset file: "+path.string());
        }

        std::string line;
        if(!std::getline(input,line))
        {
            throw std::runtime_error("dataset file is empty");
        }

        Dataset dataset;
        bool first_row=true;
        while(std::getline(input,line))
        {
            if(line.empty())
            {
                continue;
            }

            const auto fields=split_csv_line(line);
            if(fields.size()<8)
            {
                throw std::runtime_error("dataset row has too few columns");
            }

            AwgnDatasetConfig row_config;
            row_config.length=static_cast<std::size_t>(std::stoull(fields[1]));
            row_config.payload_bits=static_cast<std::size_t>(std::stoull(fields[2]));
            row_config.crc_bits=static_cast<std::size_t>(std::stoull(fields[3]));
            row_config.ebn0_db=std::stod(fields[4]);
            row_config.seed=static_cast<std::uint64_t>(std::stoull(fields[5]));
            if(fields.size()!=7+row_config.length)
            {
                throw std::runtime_error("dataset row does not match declared code length");
            }

            if(first_row)
            {
                dataset.config=row_config;
                first_row=false;
            }

            else if(row_config.length!=dataset.config.length||
            row_config.payload_bits!=dataset.config.payload_bits||
            row_config.crc_bits!=dataset.config.crc_bits||
            row_config.ebn0_db!=dataset.config.ebn0_db||
            row_config.seed!=dataset.config.seed)
            {
                throw std::runtime_error("dataset rows have inconsistent configuration metadata");
            }

            ExperimentFrame frame;
            frame.payload=string_to_bits(fields[6]);
            frame.llrs.reserve(row_config.length);
            for(std::size_t index=0;index<row_config.length;++index)
            {
                frame.llrs.push_back(std::stod(fields[7+index]));
            }

            dataset.frames.push_back(std::move(frame));
        }

        if(dataset.frames.empty())
        {
            throw std::runtime_error("dataset contains no frames");
        }

        dataset.config.frames=dataset.frames.size();
        return dataset;
    }

    void save_dataset_binary(const Dataset& dataset,const std::filesystem::path& path)
    {
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }

        std::ofstream output(path,std::ios::binary);
        if(!output)
        {
            throw std::runtime_error("failed to create binary dataset file: "+path.string());
        }

        constexpr char magic[8]={'S','C','A','D','A','T','0','1'};
        output.write(magic,8);
        const std::uint64_t length=dataset.config.length;
        const std::uint64_t payload_bits=dataset.config.payload_bits;
        const std::uint64_t crc_bits=dataset.config.crc_bits;
        const std::uint64_t frames=dataset.frames.size();
        write_binary(output,length);
        write_binary(output,payload_bits);
        write_binary(output,crc_bits);
        write_binary(output,dataset.config.ebn0_db);
        write_binary(output,frames);
        write_binary(output,dataset.config.seed);

        for(const ExperimentFrame& frame:dataset.frames)
        {
            if(frame.payload.size()!=dataset.config.payload_bits||frame.llrs.size()!=dataset.config.length)
            {
                throw std::invalid_argument("dataset frame shape does not match binary header");
            }
            for(Bit bit:frame.payload)
            {
                const std::uint8_t value=bit;
                write_binary(output,value);
            }
            output.write(
                reinterpret_cast<const char*>(frame.llrs.data()),
                static_cast<std::streamsize>(frame.llrs.size()*sizeof(double)));
            if(!output)
            {
                throw std::runtime_error("failed to write binary LLR payload");
            }
        }
    }

    Dataset load_dataset_binary(const std::filesystem::path& path)
    {
        std::ifstream input(path,std::ios::binary);
        if(!input)
        {
            throw std::runtime_error("failed to open binary dataset file: "+path.string());
        }

        char magic[8]{};
        input.read(magic,8);
        constexpr char expected[8]={'S','C','A','D','A','T','0','1'};
        if(!input||std::memcmp(magic,expected,8)!=0)
        {
            throw std::runtime_error("invalid binary dataset magic");
        }

        Dataset dataset;
        dataset.config.length=static_cast<std::size_t>(read_binary<std::uint64_t>(input));
        dataset.config.payload_bits=static_cast<std::size_t>(read_binary<std::uint64_t>(input));
        dataset.config.crc_bits=static_cast<std::size_t>(read_binary<std::uint64_t>(input));
        dataset.config.ebn0_db=read_binary<double>(input);
        dataset.config.frames=static_cast<std::size_t>(read_binary<std::uint64_t>(input));
        dataset.config.seed=read_binary<std::uint64_t>(input);
        dataset.frames.reserve(dataset.config.frames);

        for(std::size_t frame_index=0;frame_index<dataset.config.frames;++frame_index)
        {
            ExperimentFrame frame;
            frame.payload.resize(dataset.config.payload_bits);
            for(Bit& bit:frame.payload)
            {
                bit=static_cast<Bit>(read_binary<std::uint8_t>(input));
                if(bit>1)
                {
                    throw std::runtime_error("invalid payload bit in binary dataset");
                }
            }
            frame.llrs.resize(dataset.config.length);
            input.read(
                reinterpret_cast<char*>(frame.llrs.data()),
                static_cast<std::streamsize>(frame.llrs.size()*sizeof(double)));
            if(!input)
            {
                throw std::runtime_error("unexpected end of binary LLR payload");
            }
            dataset.frames.push_back(std::move(frame));
        }
        return dataset;
    }

    UnifiedDecodeResult decode_variant(
    DecoderVariant variant,
    const std::vector<double>& llrs,
    const std::vector<bool>& frozen,
    std::size_t list_size,
    const PolarAlgorithmIR& portfolio_algorithm)
    {
        if(variant!=DecoderVariant::SC&&variant!=DecoderVariant::SCL)
        {
            const ParentPath probe{0,0.0,llrs};
            if(rate1_metric_tie_risk({probe}))
            {
                auto result=scl_decode(llrs,frozen,list_size);
                return {std::move(result.paths),result.stats,{},1};
            }
        }

        switch(variant)
        {
            case DecoderVariant::SC:
            {
                auto result=scl_decode(llrs,frozen,1);
                return {std::move(result.paths),result.stats,
                    {
                    }

                    ,0};
            }

            case DecoderVariant::SCL:
            {
                auto result=scl_decode(llrs,frozen,list_size);
                return {std::move(result.paths),result.stats,
                    {
                    }

                    ,0};
            }

            case DecoderVariant::SSCLRate01Exact:
            {
                auto result=portfolio_scl_decode(llrs,frozen,list_size,fixed_algorithm(Rate1Engine::LazyExact));
                return {
                    std::move(result.paths),
                    result.stats.decoder,
                    std::move(result.stats.engine_counts),
                    result.stats.tie_fallbacks
                };
            }

            case DecoderVariant::FastSSCLRate1Bounded:
            {
                auto result=portfolio_scl_decode(llrs,frozen,list_size,fixed_algorithm(Rate1Engine::FastSSCLBounded));
                return {
                    std::move(result.paths),
                    result.stats.decoder,
                    std::move(result.stats.engine_counts),
                    result.stats.tie_fallbacks
                };
            }

            case DecoderVariant::CertifiedFrontier:
            {
                auto result=portfolio_scl_decode(
                    llrs,
                    frozen,
                    list_size,
                    fixed_algorithm(Rate1Engine::CertifiedFrontier));
                return {
                    std::move(result.paths),
                    result.stats.decoder,
                    std::move(result.stats.engine_counts),
                    result.stats.tie_fallbacks
                };
            }

            case DecoderVariant::PortfolioRule:
            {
                auto result=portfolio_scl_decode(llrs,frozen,list_size,portfolio_algorithm);
                return {
                    std::move(result.paths),
                    result.stats.decoder,
                    std::move(result.stats.engine_counts),
                    result.stats.tie_fallbacks
                };
            }

            case DecoderVariant::SearchedPortfolio:
            {
                const SelectorGenome& genome=portfolio_algorithm.selector;
                std::optional<Rate1Engine> fixed_engine;
                if(genome.list_size_threshold>0&&list_size<=genome.list_size_threshold)
                {
                    fixed_engine=genome.low_list_engine;
                }
                else if(genome.small_engine==genome.easy_engine&&genome.easy_engine==genome.hard_engine)
                {
                    fixed_engine=genome.small_engine;
                }

                const bool has_special_nodes=genome.enable_repetition||genome.enable_spc;
                if(!has_special_nodes&&fixed_engine)
                {
                    auto result=portfolio_scl_decode(llrs,frozen,list_size,fixed_algorithm(*fixed_engine));
                    return {
                        std::move(result.paths),
                        result.stats.decoder,
                        std::move(result.stats.engine_counts),
                        result.stats.tie_fallbacks
                    };
                }

                auto result=portfolio_scl_decode(llrs,frozen,list_size,portfolio_algorithm);
                return {
                    std::move(result.paths),
                    result.stats.decoder,
                    std::move(result.stats.engine_counts),
                    result.stats.tie_fallbacks
                };
            }
        }

        throw std::logic_error("unreachable decoder variant");
    }

    bool paths_equivalent(
    const std::vector<DecoderPath>& lhs,
    const std::vector<DecoderPath>& rhs,
    double tolerance)
    {
        if(lhs.size()!=rhs.size())
        {
            return false;
        }

        for(std::size_t index=0;index<lhs.size();++index)
        {
            if(lhs[index].u!=rhs[index].u)
            {
                return false;
            }

            const double scale=std::max({1.0,std::abs(lhs[index].metric),std::abs(rhs[index].metric)});
            if(std::abs(lhs[index].metric-rhs[index].metric)>tolerance*scale)
            {
                return false;
            }
        }

        return true;
    }

    double ExperimentMetrics::fer() const
    {
        return frames==0?0.0:static_cast<double>(frame_errors)/static_cast<double>(frames);
    }

    double ExperimentMetrics::ber() const
    {
        const std::size_t total_bits=frames*payload_bits;
        return total_bits==0?0.0:static_cast<double>(bit_errors)/static_cast<double>(total_bits);
    }

    double ExperimentMetrics::mismatch_rate() const
    {
        return compared_frames==0?0.0:static_cast<double>(list_mismatches)/static_cast<double>(compared_frames);
    }

    ExperimentMetrics run_decoder_experiment(
    const Dataset& dataset,
    DecoderVariant variant,
    const ExperimentRunConfig& run_config,
    const PolarAlgorithmIR& portfolio_algorithm)
    {
        if(run_config.list_size==0)
        {
            throw std::invalid_argument("list_size must be positive");
        }

        if(run_config.timing_repeats==0)
        {
            throw std::invalid_argument("timing_repeats must be positive");
        }

        const auto frozen=pw_frozen_pattern(dataset.config.length,dataset.config.information_bits());
        ExperimentMetrics metrics;
        metrics.variant=variant;
        metrics.length=dataset.config.length;
        metrics.information_bits=dataset.config.information_bits();
        metrics.payload_bits=dataset.config.payload_bits;
        metrics.crc_bits=dataset.config.crc_bits;
        metrics.list_size=variant==DecoderVariant::SC?1:run_config.list_size;
        metrics.ebn0_db=dataset.config.ebn0_db;
        metrics.frames=dataset.frames.size();

        const std::size_t warmup_count=std::min(run_config.warmup_frames,dataset.frames.size());
        for(std::size_t index=0;index<warmup_count;++index)
        {
            const auto warmup=decode_variant(
            variant,
            dataset.frames[index].llrs,
            frozen,
            run_config.list_size,
            portfolio_algorithm);
            if(warmup.paths.empty())
            {
                throw std::runtime_error("decoder returned an empty path list during warmup");
            }
        }

        std::vector<double> latencies;
        latencies.reserve(dataset.frames.size()*run_config.timing_repeats);
        double leaf_sum=0.0;
        double split_sum=0.0;
        double prune_sum=0.0;
        double rate0_sum=0.0;
        double rate1_sum=0.0;
        double repetition_sum=0.0;
        double spc_sum=0.0;
        double special_candidate_sum=0.0;
        double candidate_sum=0.0;
        double expansion_sum=0.0;
        double fallback_sum=0.0;
        double tie_fallback_sum=0.0;
        double exhaustive_sum=0.0;
        double lazy_sum=0.0;
        double fast_sscl_sum=0.0;
        double frontier_sum=0.0;

        for(const auto& frame:dataset.frames)
        {
            UnifiedDecodeResult result;
            for(std::size_t repeat=0;repeat<run_config.timing_repeats;++repeat)
            {
                const auto start=std::chrono::steady_clock::now();
                auto timed=decode_variant(
                variant,
                frame.llrs,
                frozen,
                run_config.list_size,
                portfolio_algorithm);
                const auto stop=std::chrono::steady_clock::now();
                latencies.push_back(std::chrono::duration<double,std::micro>(stop-start).count());
                if(repeat==0)
                {
                    result=std::move(timed);
                }
            }

            if(result.paths.empty())
            {
                throw std::runtime_error("decoder returned an empty path list");
            }

            DecoderPath selected;
            if(variant==DecoderVariant::SC)
            {
                selected=result.paths.front();
            }

            else
            {
                selected=select_crc_path(result.paths,frozen,dataset.config.crc_bits);
            }

            const auto information=extract_information(selected.u,frozen);
            const std::size_t frame_bit_errors=count_payload_errors(frame.payload,information);
            metrics.bit_errors+=frame_bit_errors;
            metrics.frame_errors+=frame_bit_errors>0?1U:0U;

            if(run_config.compare_with_scl&&variant!=DecoderVariant::SC&&variant!=DecoderVariant::SCL)
            {
                const auto reference=scl_decode(frame.llrs,frozen,run_config.list_size);
                ++metrics.compared_frames;
                metrics.list_mismatches+=paths_equivalent(reference.paths,result.paths)?0U:1U;
            }

            leaf_sum+=static_cast<double>(result.stats.leaf_steps);
            split_sum+=static_cast<double>(result.stats.path_splits);
            prune_sum+=static_cast<double>(result.stats.sorts);
            rate0_sum+=static_cast<double>(result.stats.rate0_nodes);
            rate1_sum+=static_cast<double>(result.stats.rate1_nodes);
            repetition_sum+=static_cast<double>(result.stats.repetition_nodes);
            spc_sum+=static_cast<double>(result.stats.spc_nodes);
            special_candidate_sum+=static_cast<double>(result.stats.special_candidates);
            candidate_sum+=static_cast<double>(result.stats.rate1_candidates);
            expansion_sum+=static_cast<double>(result.stats.rate1_frontier_expansions);
            fallback_sum+=static_cast<double>(result.stats.rate1_fallbacks);
            tie_fallback_sum+=static_cast<double>(result.tie_fallbacks);
            const auto count_engine=[&](Rate1Engine engine)
            {
                const auto iterator=result.engine_counts.find(engine);
                return iterator==result.engine_counts.end()?0.0:static_cast<double>(iterator->second);
            };
            exhaustive_sum+=count_engine(Rate1Engine::Exhaustive);
            lazy_sum+=count_engine(Rate1Engine::LazyExact);
            fast_sscl_sum+=count_engine(Rate1Engine::FastSSCLBounded);
            frontier_sum+=count_engine(Rate1Engine::CertifiedFrontier);
        }

        const double frame_count=static_cast<double>(dataset.frames.size());
        metrics.mean_leaf_steps=leaf_sum/frame_count;
        metrics.mean_path_splits=split_sum/frame_count;
        metrics.mean_prune_events=prune_sum/frame_count;
        metrics.mean_rate0_nodes=rate0_sum/frame_count;
        metrics.mean_rate1_nodes=rate1_sum/frame_count;
        metrics.mean_repetition_nodes=repetition_sum/frame_count;
        metrics.mean_spc_nodes=spc_sum/frame_count;
        metrics.mean_special_candidates=special_candidate_sum/frame_count;
        metrics.mean_rate1_candidates=candidate_sum/frame_count;
        metrics.mean_frontier_expansions=expansion_sum/frame_count;
        metrics.mean_rate1_fallbacks=fallback_sum/frame_count;
        metrics.mean_tie_fallbacks=tie_fallback_sum/frame_count;
        metrics.mean_exhaustive_nodes=exhaustive_sum/frame_count;
        metrics.mean_lazy_exact_nodes=lazy_sum/frame_count;
        metrics.mean_fast_sscl_nodes=fast_sscl_sum/frame_count;
        metrics.mean_certified_frontier_nodes=frontier_sum/frame_count;
        metrics.mean_us=std::accumulate(latencies.begin(),latencies.end(),0.0)/static_cast<double>(latencies.size());
        metrics.p50_us=percentile(latencies,0.50);
        metrics.p95_us=percentile(latencies,0.95);
        metrics.p99_us=percentile(latencies,0.99);
        metrics.max_us=*std::max_element(latencies.begin(),latencies.end());
        metrics.throughput_frames_per_second=metrics.mean_us>0.0?1.0e6/metrics.mean_us:0.0;
        return metrics;
    }

    void write_metrics_csv(
    const std::vector<ExperimentMetrics>& metrics,
    const std::filesystem::path& path)
    {
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }

        std::ofstream output(path);
        if(!output)
        {
            throw std::runtime_error("failed to create metrics CSV: "+path.string());
        }

        output<<"decoder,N,K,payload_bits,crc_bits,L,EbN0_dB,frames,frame_errors,bit_errors,FER,BER,compared_frames,list_mismatches,mismatch_rate,mean_us,p50_us,p95_us,p99_us,max_us,throughput_fps,mean_leaf_steps,mean_path_splits,mean_prune_events,mean_rate0_nodes,mean_rate1_nodes,mean_repetition_nodes,mean_spc_nodes,mean_special_candidates,mean_rate1_candidates,mean_frontier_expansions,mean_rate1_fallbacks,mean_tie_fallbacks,mean_exhaustive_nodes,mean_lazy_exact_nodes,mean_fast_sscl_nodes,mean_certified_frontier_nodes\n";
        output<<std::setprecision(12);
        for(const auto& row:metrics)
        {
            output<<to_string(row.variant)<<','
            <<row.length<<','
            <<row.information_bits<<','
            <<row.payload_bits<<','
            <<row.crc_bits<<','
            <<row.list_size<<','
            <<row.ebn0_db<<','
            <<row.frames<<','
            <<row.frame_errors<<','
            <<row.bit_errors<<','
            <<row.fer()<<','
            <<row.ber()<<','
            <<row.compared_frames<<','
            <<row.list_mismatches<<','
            <<row.mismatch_rate()<<','
            <<row.mean_us<<','
            <<row.p50_us<<','
            <<row.p95_us<<','
            <<row.p99_us<<','
            <<row.max_us<<','
            <<row.throughput_frames_per_second<<','
            <<row.mean_leaf_steps<<','
            <<row.mean_path_splits<<','
            <<row.mean_prune_events<<','
            <<row.mean_rate0_nodes<<','
            <<row.mean_rate1_nodes<<','
            <<row.mean_repetition_nodes<<','
            <<row.mean_spc_nodes<<','
            <<row.mean_special_candidates<<','
            <<row.mean_rate1_candidates<<','
            <<row.mean_frontier_expansions<<','
            <<row.mean_rate1_fallbacks<<','
            <<row.mean_tie_fallbacks<<','
            <<row.mean_exhaustive_nodes<<','
            <<row.mean_lazy_exact_nodes<<','
            <<row.mean_fast_sscl_nodes<<','
            <<row.mean_certified_frontier_nodes<<'\n';
        }
    }
}
// namespace polar_scl
