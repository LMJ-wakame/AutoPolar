#include "polar_scl/polar_scl.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace polar_scl;

namespace
{
    struct SuiteOptions
    {
        std::string profile{"quick"};
        std::size_t frames{0};
        std::uint64_t seed_offset{300000000ULL};
        std::filesystem::path output_directory{"results/paper"};
        std::optional<std::filesystem::path> selector_path;
        bool save_datasets{true};
    };

    SuiteOptions parse_options(int argc,char** argv)
    {
        SuiteOptions options;
        for(int index=1;index<argc;++index)
        {
            const std::string argument=argv[index];
            if(argument=="--profile"&&index+1<argc)
            {
                options.profile=argv[++index];
            }
            else if(argument=="--frames"&&index+1<argc)
            {
                options.frames=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--seed-offset"&&index+1<argc)
            {
                options.seed_offset=static_cast<std::uint64_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--out"&&index+1<argc)
            {
                options.output_directory=argv[++index];
            }
            else if(argument=="--selector"&&index+1<argc)
            {
                options.selector_path=std::filesystem::path(argv[++index]);
            }
            else if(argument=="--no-save-datasets")
            {
                options.save_datasets=false;
            }
            else
            {
                throw std::invalid_argument("unknown argument: "+argument);
            }
        }

        if(options.profile!="quick"&&options.profile!="paper")
        {
            throw std::invalid_argument("profile must be quick or paper");
        }

        return options;
    }

    std::vector<std::size_t> profile_lengths(const SuiteOptions& options)
    {
        return options.profile=="quick"?std::vector<std::size_t>{64,128}
            :std::vector<std::size_t>{128,256,512};
    }

    std::vector<double> profile_rates()
    {
        return {0.5,0.75};
    }

    std::vector<double> profile_ebn0(const SuiteOptions& options)
    {
        return options.profile=="quick"?std::vector<double>{0.5,1.5,2.5}
            :std::vector<double>{0.0,0.5,1.0,1.5,2.0,2.5,3.0};
    }

    std::vector<std::size_t> profile_list_sizes(const SuiteOptions& options)
    {
        return options.profile=="quick"?std::vector<std::size_t>{4,8}
            :std::vector<std::size_t>{4,8,16};
    }

    std::size_t profile_frames(const SuiteOptions& options)
    {
        if(options.frames!=0)
        {
            return options.frames;
        }

        return options.profile=="quick"?100U:1000U;
    }

    std::uint64_t dataset_seed(
        const SuiteOptions& options,
        std::size_t length,
        std::size_t information_bits,
        double ebn0)
    {
        return 202608160000ULL+options.seed_offset
            +static_cast<std::uint64_t>(length)*10000ULL
            +static_cast<std::uint64_t>(information_bits)*100ULL
            +static_cast<std::uint64_t>(std::llround(ebn0*10.0));
    }

    std::string dataset_filename(const AwgnDatasetConfig& config)
    {
        std::ostringstream stream;
        stream<<"N"<<config.length
            <<"_K"<<config.information_bits()
            <<"_EbN0_"<<std::fixed<<std::setprecision(1)<<config.ebn0_db
            <<"dB.csv";
        return stream.str();
    }

    std::pair<double,double> wilson95(std::size_t errors,std::size_t total)
    {
        if(total==0)
        {
            return {0.0,0.0};
        }

        constexpr double z=1.959963984540054;
        const double n=static_cast<double>(total);
        const double p=static_cast<double>(errors)/n;
        const double z2=z*z;
        const double denominator=1.0+z2/n;
        const double center=(p+z2/(2.0*n))/denominator;
        const double half=z*std::sqrt((p*(1.0-p)+z2/(4.0*n))/n)/denominator;
        return {std::max(0.0,center-half),std::min(1.0,center+half)};
    }

    void write_fer_ber(const std::vector<ExperimentMetrics>& rows,const std::filesystem::path& path)
    {
        std::ofstream output(path);
        output<<"decoder,N,K,R,L,EbN0_dB,frames,frame_errors,bit_errors,FER,FER95_low,FER95_high,BER,list_mismatches,mismatch_rate\n";
        output<<std::setprecision(12);
        for(const auto& row:rows)
        {
            const double rate=static_cast<double>(row.information_bits)/static_cast<double>(row.length);
            const auto [fer_low,fer_high]=wilson95(row.frame_errors,row.frames);
            output<<to_string(row.variant)<<','
                <<row.length<<','
                <<row.information_bits<<','
                <<rate<<','
                <<row.list_size<<','
                <<row.ebn0_db<<','
                <<row.frames<<','
                <<row.frame_errors<<','
                <<row.bit_errors<<','
                <<row.fer()<<','
                <<fer_low<<','
                <<fer_high<<','
                <<row.ber()<<','
                <<row.list_mismatches<<','
                <<row.mismatch_rate()<<'\n';
        }
    }

    void write_latency(const std::vector<ExperimentMetrics>& rows,const std::filesystem::path& path)
    {
        std::ofstream output(path);
        output<<"decoder,N,K,L,EbN0_dB,frames,mean_us,p50_us,p95_us,p99_us,max_us,throughput_fps\n";
        output<<std::setprecision(12);
        for(const auto& row:rows)
        {
            output<<to_string(row.variant)<<','
                <<row.length<<','
                <<row.information_bits<<','
                <<row.list_size<<','
                <<row.ebn0_db<<','
                <<row.frames<<','
                <<row.mean_us<<','
                <<row.p50_us<<','
                <<row.p95_us<<','
                <<row.p99_us<<','
                <<row.max_us<<','
                <<row.throughput_frames_per_second<<'\n';
        }
    }

    void write_complexity(const std::vector<ExperimentMetrics>& rows,const std::filesystem::path& path)
    {
        std::ofstream output(path);
        output<<"decoder,N,K,L,EbN0_dB,mean_leaf_steps,mean_path_splits,mean_prune_events,mean_rate0_nodes,mean_rate1_nodes,mean_repetition_nodes,mean_spc_nodes,mean_special_candidates,mean_rate1_candidates,mean_frontier_expansions,mean_rate1_fallbacks,mean_tie_fallbacks,mean_exhaustive_nodes,mean_lazy_exact_nodes,mean_fast_sscl_nodes,mean_certified_frontier_nodes\n";
        output<<std::setprecision(12);
        for(const auto& row:rows)
        {
            output<<to_string(row.variant)<<','
                <<row.length<<','
                <<row.information_bits<<','
                <<row.list_size<<','
                <<row.ebn0_db<<','
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

    using ResultKey=std::tuple<std::size_t,std::size_t,std::size_t,long long>;

    ResultKey make_key(const ExperimentMetrics& row)
    {
        return {
            row.length,
            row.information_bits,
            row.list_size,
            static_cast<long long>(std::llround(row.ebn0_db*1000.0))};
    }

    void write_comparison(const std::vector<ExperimentMetrics>& rows,const std::filesystem::path& path)
    {
        std::map<ResultKey,const ExperimentMetrics*> scl_rows;
        for(const auto& row:rows)
        {
            if(row.variant==DecoderVariant::SCL)
            {
                scl_rows[make_key(row)]=&row;
            }
        }

        std::ofstream output(path);
        output<<"decoder,N,K,L,EbN0_dB,FER,FER_delta_vs_SCL,mean_us,speedup_vs_SCL,latency_reduction_percent,p99_us,mismatch_rate\n";
        output<<std::setprecision(12);
        for(const auto& row:rows)
        {
            const auto found=scl_rows.find(make_key(row));
            if(found==scl_rows.end())
            {
                continue;
            }

            const auto& reference=*found->second;
            const double speedup=row.mean_us>0.0?reference.mean_us/row.mean_us:0.0;
            const double reduction=reference.mean_us>0.0?100.0*(1.0-row.mean_us/reference.mean_us):0.0;
            output<<to_string(row.variant)<<','
                <<row.length<<','
                <<row.information_bits<<','
                <<row.list_size<<','
                <<row.ebn0_db<<','
                <<row.fer()<<','
                <<row.fer()-reference.fer()<<','
                <<row.mean_us<<','
                <<speedup<<','
                <<reduction<<','
                <<row.p99_us<<','
                <<row.mismatch_rate()<<'\n';
        }
    }

    void write_ablation(const std::vector<ExperimentMetrics>& rows,const std::filesystem::path& path)
    {
        std::ofstream output(path);
        output<<"variant,N,K,L,EbN0_dB,mean_us,p99_us,FER,mean_candidates,frontier_expansions,fallbacks,rep_nodes,spc_nodes,special_candidates,lazy_nodes,fast_sscl_nodes,frontier_nodes\n";
        output<<std::setprecision(12);
        for(const auto& row:rows)
        {
            if(row.variant==DecoderVariant::SC||row.variant==DecoderVariant::SCL)
            {
                continue;
            }

            output<<to_string(row.variant)<<','
                <<row.length<<','
                <<row.information_bits<<','
                <<row.list_size<<','
                <<row.ebn0_db<<','
                <<row.mean_us<<','
                <<row.p99_us<<','
                <<row.fer()<<','
                <<row.mean_rate1_candidates<<','
                <<row.mean_frontier_expansions<<','
                <<row.mean_rate1_fallbacks<<','
                <<row.mean_repetition_nodes<<','
                <<row.mean_spc_nodes<<','
                <<row.mean_special_candidates<<','
                <<row.mean_lazy_exact_nodes<<','
                <<row.mean_fast_sscl_nodes<<','
                <<row.mean_certified_frontier_nodes<<'\n';
        }
    }


    void write_special_node_ablation(
        const std::vector<std::pair<ExperimentMetrics,ExperimentMetrics>>& pairs,
        const std::filesystem::path& path)
    {
        std::ofstream output(path);
        output<<"N,K,L,EbN0_dB,FER_no_special,FER_v2,mean_us_no_special,mean_us_v2,speedup_from_special_nodes,latency_reduction_percent,rep_nodes,spc_nodes,special_candidates,mismatch_v2\n";
        output<<std::setprecision(12);
        for(const auto& [without_special,with_special]:pairs)
        {
            const double speedup=with_special.mean_us>0.0?without_special.mean_us/with_special.mean_us:0.0;
            const double reduction=without_special.mean_us>0.0?
                100.0*(1.0-with_special.mean_us/without_special.mean_us):0.0;
            output<<with_special.length<<','
                <<with_special.information_bits<<','
                <<with_special.list_size<<','
                <<with_special.ebn0_db<<','
                <<without_special.fer()<<','
                <<with_special.fer()<<','
                <<without_special.mean_us<<','
                <<with_special.mean_us<<','
                <<speedup<<','
                <<reduction<<','
                <<with_special.mean_repetition_nodes<<','
                <<with_special.mean_spc_nodes<<','
                <<with_special.mean_special_candidates<<','
                <<with_special.mismatch_rate()<<'\n';
        }
    }

    void write_aggregate_summary(const std::vector<ExperimentMetrics>& rows,const std::filesystem::path& path)
    {
        std::map<ResultKey,const ExperimentMetrics*> scl_rows;
        for(const auto& row:rows)
        {
            if(row.variant==DecoderVariant::SCL)
            {
                scl_rows[make_key(row)]=&row;
            }
        }

        struct Aggregate
        {
            std::size_t configurations{0};
            double speedup_sum{0.0};
            double reduction_sum{0.0};
            double fer_delta_sum{0.0};
            std::size_t mismatches{0};
            std::size_t compared_frames{0};
            std::vector<double> speedups;
        };

        std::map<DecoderVariant,Aggregate> aggregates;
        for(const auto& row:rows)
        {
            if(row.variant==DecoderVariant::SC||row.variant==DecoderVariant::SCL)
            {
                continue;
            }

            const auto reference_iterator=scl_rows.find(make_key(row));
            if(reference_iterator==scl_rows.end())
            {
                continue;
            }

            const ExperimentMetrics& reference=*reference_iterator->second;
            const double speedup=row.mean_us>0.0?reference.mean_us/row.mean_us:0.0;
            const double reduction=reference.mean_us>0.0?100.0*(1.0-row.mean_us/reference.mean_us):0.0;
            Aggregate& aggregate=aggregates[row.variant];
            ++aggregate.configurations;
            aggregate.speedup_sum+=speedup;
            aggregate.reduction_sum+=reduction;
            aggregate.fer_delta_sum+=row.fer()-reference.fer();
            aggregate.mismatches+=row.list_mismatches;
            aggregate.compared_frames+=row.compared_frames;
            aggregate.speedups.push_back(speedup);
        }

        std::ofstream output(path);
        output<<"decoder,configurations,mean_speedup_vs_SCL,median_speedup_vs_SCL,mean_latency_reduction_percent,mean_FER_delta,total_list_mismatches,total_compared_frames\n";
        output<<std::setprecision(12);
        for(auto& [variant,aggregate]:aggregates)
        {
            std::sort(aggregate.speedups.begin(),aggregate.speedups.end());
            double median=0.0;
            if(!aggregate.speedups.empty())
            {
                const std::size_t middle=aggregate.speedups.size()/2;
                median=aggregate.speedups.size()%2==0
                    ?(aggregate.speedups[middle-1]+aggregate.speedups[middle])/2.0
                    :aggregate.speedups[middle];
            }

            const double count=static_cast<double>(aggregate.configurations);
            output<<to_string(variant)<<','
                <<aggregate.configurations<<','
                <<(count>0.0?aggregate.speedup_sum/count:0.0)<<','
                <<median<<','
                <<(count>0.0?aggregate.reduction_sum/count:0.0)<<','
                <<(count>0.0?aggregate.fer_delta_sum/count:0.0)<<','
                <<aggregate.mismatches<<','
                <<aggregate.compared_frames<<'\n';
        }
    }

    void write_summary(
        const std::vector<ExperimentMetrics>& rows,
        const SuiteOptions& options,
        bool searched_enabled,
        const std::filesystem::path& path)
    {
        std::ofstream output(path);
        output<<"# Polar Fast-SCL paper experiment summary\n\n";
        output<<"Profile: `"<<options.profile<<"`  \n";
        output<<"Rows: "<<rows.size()<<"  \n";
        output<<"Test seed offset: "<<options.seed_offset<<"  \n";
        output<<"Construction: polarization-weight (PW), not full 5G NR construction.  \n";
        output<<"Channel: BPSK over AWGN; noise variance uses Eb/N0 and R=K/N.  \n";
        output<<"All decoders in a configuration consume the same generated LLR frames.  \n\n";
        output<<"## Methods\n\n";
        output<<"- SC: bit-wise successive cancellation (`L=1`).\n";
        output<<"- SCL: bit-wise list decoder used as the exact reference.\n";
        output<<"- SSCL-Rate01-Exact: maximal Rate-0/Rate-1 schedule with exact lazy Rate-1 Top-L.\n";
        output<<"- FastSSCL-Rate1-Bounded: software realization of the Rate-1 list-size reliability bound.\n";
        output<<"- Certified-Frontier: adaptive frontier with certificate and exact fallback.\n";
        output<<"- Portfolio-Rule: hand-designed Algorithm-IR rule selector.\n";
        if(searched_enabled)
        {
            output<<"- SCA-Polar: searched SelectorGenome selected on disjoint train/validation splits, evaluated only here on the independent test split.\n";
        }
        output<<"\nThe quick profile is a pipeline smoke test. Use much larger frame counts for final low-FER claims.\n";
    }
}

int main(int argc,char** argv)
{
    try
    {
        const SuiteOptions options=parse_options(argc,argv);
        const std::size_t frames=profile_frames(options);
        std::filesystem::create_directories(options.output_directory);
        const auto dataset_directory=options.output_directory/"datasets";
        if(options.save_datasets)
        {
            std::filesystem::create_directories(dataset_directory);
        }

        PolarAlgorithmIR searched_algorithm;
        const bool searched_enabled=options.selector_path.has_value();
        if(searched_enabled)
        {
            searched_algorithm.selector_mode=SelectorMode::Rule;
            searched_algorithm.selector=load_selector_genome(*options.selector_path);
            validate_algorithm_ir(searched_algorithm);
        }

        std::cout<<"Polar paper suite profile="<<options.profile
            <<" frames/config="<<frames
            <<" output="<<options.output_directory.string()
            <<" searched="<<(searched_enabled?"yes":"no")<<'\n';

        std::vector<ExperimentMetrics> rows;
        std::vector<std::pair<ExperimentMetrics,ExperimentMetrics>> special_ablation_pairs;
        for(std::size_t length:profile_lengths(options))
        {
            for(double rate:profile_rates())
            {
                const std::size_t information_bits=static_cast<std::size_t>(
                    std::llround(rate*static_cast<double>(length)));
                constexpr std::size_t crc_bits=8;
                if(information_bits<=crc_bits)
                {
                    continue;
                }

                const std::size_t payload_bits=information_bits-crc_bits;
                for(double ebn0:profile_ebn0(options))
                {
                    AwgnDatasetConfig dataset_config;
                    dataset_config.length=length;
                    dataset_config.payload_bits=payload_bits;
                    dataset_config.crc_bits=crc_bits;
                    dataset_config.ebn0_db=ebn0;
                    dataset_config.frames=frames;
                    dataset_config.seed=dataset_seed(options,length,information_bits,ebn0);
                    const Dataset dataset=generate_awgn_dataset(dataset_config);
                    if(options.save_datasets)
                    {
                        save_dataset_csv(dataset,dataset_directory/dataset_filename(dataset_config));
                    }

                    for(std::size_t list_size:profile_list_sizes(options))
                    {
                        ExperimentRunConfig run_config;
                        run_config.list_size=list_size;
                        run_config.compare_with_scl=true;
                        run_config.warmup_frames=std::min<std::size_t>(8,frames);
                        run_config.timing_repeats=1;
                        for(DecoderVariant variant:paper_decoder_variants(searched_enabled))
                        {
                            const PolarAlgorithmIR algorithm=variant==DecoderVariant::SearchedPortfolio
                                ?searched_algorithm:PolarAlgorithmIR{};
                            ExperimentMetrics metrics=run_decoder_experiment(
                                dataset,variant,run_config,algorithm);
                            std::cout<<"[RESULT] "<<std::setw(24)<<std::left<<to_string(variant)
                                <<std::right
                                <<" N="<<length
                                <<" K="<<information_bits
                                <<" L="<<metrics.list_size
                                <<" EbN0="<<std::fixed<<std::setprecision(1)<<ebn0
                                <<" FER="<<std::setprecision(4)<<metrics.fer()
                                <<" mean_us="<<std::setprecision(2)<<metrics.mean_us
                                <<" mismatch="<<metrics.list_mismatches<<'\n';
                            if(variant==DecoderVariant::SearchedPortfolio)
                            {
                                PolarAlgorithmIR no_special=searched_algorithm;
                                no_special.selector.enable_repetition=false;
                                no_special.selector.enable_spc=false;
                                ExperimentMetrics without_special=run_decoder_experiment(
                                    dataset,variant,run_config,no_special);
                                special_ablation_pairs.push_back({without_special,metrics});
                            }
                            rows.push_back(std::move(metrics));
                        }
                    }
                }
            }
        }

        write_metrics_csv(rows,options.output_directory/"all_metrics.csv");
        write_fer_ber(rows,options.output_directory/"fer_ber.csv");
        write_latency(rows,options.output_directory/"latency.csv");
        write_complexity(rows,options.output_directory/"complexity.csv");
        write_comparison(rows,options.output_directory/"comparison_vs_scl.csv");
        write_ablation(rows,options.output_directory/"ablation.csv");
        if(!special_ablation_pairs.empty())
        {
            write_special_node_ablation(
                special_ablation_pairs,options.output_directory/"special_node_ablation.csv");
        }
        write_aggregate_summary(rows,options.output_directory/"aggregate_vs_scl.csv");
        write_summary(rows,options,searched_enabled,options.output_directory/"SUMMARY.md");
        std::cout<<"Generated paper CSV files in "<<options.output_directory.string()<<'\n';
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"[ERROR] "<<error.what()<<'\n';
        return 1;
    }
}
