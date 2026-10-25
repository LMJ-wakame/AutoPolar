#include "polar_scl/polar_scl.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace polar_scl;

namespace
{
    struct CorrectnessRow
    {
        std::string decoder;
        std::string llr_mode;
        std::size_t length{};
        std::size_t information_bits{};
        std::size_t list_size{};
        std::size_t trials{};
        std::size_t mismatches{};
    };

    void write_rows(const std::vector<CorrectnessRow>& rows,const std::filesystem::path& path)
    {
        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }

        std::ofstream output(path);
        if(!output)
        {
            throw std::runtime_error("cannot create correctness CSV");
        }

        output<<"decoder,llr_mode,N,K,L,trials,mismatches,mismatch_rate\n";
        output<<std::setprecision(12);
        for(const auto& row:rows)
        {
            const double rate=row.trials==0?0.0:static_cast<double>(row.mismatches)/static_cast<double>(row.trials);
            output<<row.decoder<<','
                <<row.llr_mode<<','
                <<row.length<<','
                <<row.information_bits<<','
                <<row.list_size<<','
                <<row.trials<<','
                <<row.mismatches<<','
                <<rate<<'\n';
        }
    }
}

int main(int argc,char** argv)
{
    try
    {
        std::size_t trials=80;
        std::filesystem::path output="results/paper/correctness_matrix.csv";
        std::optional<std::filesystem::path> selector_path;
        for(int index=1;index<argc;++index)
        {
            const std::string argument=argv[index];
            if(argument=="--trials"&&index+1<argc)
            {
                trials=static_cast<std::size_t>(std::stoull(argv[++index]));
            }
            else if(argument=="--out"&&index+1<argc)
            {
                output=argv[++index];
            }
            else if(argument=="--selector"&&index+1<argc)
            {
                selector_path=std::filesystem::path(argv[++index]);
            }
            else
            {
                throw std::invalid_argument("unknown argument: "+argument);
            }
        }

        PolarAlgorithmIR searched_algorithm;
        if(selector_path)
        {
            searched_algorithm.selector_mode=SelectorMode::Rule;
            searched_algorithm.selector=load_selector_genome(*selector_path);
        }

        std::vector<DecoderVariant> variants{
            DecoderVariant::SSCLRate01Exact,
            DecoderVariant::FastSSCLRate1Bounded,
            DecoderVariant::CertifiedFrontier,
            DecoderVariant::PortfolioRule};
        if(selector_path)
        {
            variants.push_back(DecoderVariant::SearchedPortfolio);
        }

        std::mt19937_64 rng(20260816);
        std::normal_distribution<double> continuous(0.0,3.0);
        std::uniform_int_distribution<int> quantized(-12,12);
        std::vector<CorrectnessRow> rows;

        for(const std::string mode:{std::string("continuous"),std::string("quantized_integer")})
        {
            for(std::size_t length:{32U,64U,128U})
            {
                for(double rate:{0.5,0.75})
                {
                    const std::size_t information_bits=static_cast<std::size_t>(
                        std::llround(rate*static_cast<double>(length)));
                    const auto frozen=pw_frozen_pattern(length,information_bits);
                    for(std::size_t list_size:{2U,4U,8U})
                    {
                        std::vector<CorrectnessRow> case_rows;
                        for(DecoderVariant variant:variants)
                        {
                            case_rows.push_back({
                                to_string(variant),mode,length,information_bits,list_size,trials,0});
                        }

                        for(std::size_t trial=0;trial<trials;++trial)
                        {
                            std::vector<double> llrs(length);
                            for(double& llr:llrs)
                            {
                                llr=mode=="continuous"?continuous(rng):static_cast<double>(quantized(rng));
                            }

                            const auto reference=scl_decode(llrs,frozen,list_size);
                            for(std::size_t variant_index=0;variant_index<variants.size();++variant_index)
                            {
                                const DecoderVariant variant=variants[variant_index];
                                const PolarAlgorithmIR algorithm=variant==DecoderVariant::SearchedPortfolio
                                    ?searched_algorithm:PolarAlgorithmIR{};
                                const auto actual=decode_variant(
                                    variant,llrs,frozen,list_size,algorithm);
                                case_rows[variant_index].mismatches+=paths_equivalent(
                                    reference.paths,
                                    actual.paths,
                                    mode=="continuous"?1e-10:0.0)?0U:1U;
                            }
                        }

                        for(const CorrectnessRow& row:case_rows)
                        {
                            rows.push_back(row);
                            std::cout<<"[CHECK] "<<row.decoder
                                <<" mode="<<row.llr_mode
                                <<" N="<<length
                                <<" K="<<information_bits
                                <<" L="<<list_size
                                <<" mismatch="<<row.mismatches<<'/'<<row.trials<<'\n';
                        }
                    }
                }
            }
        }

        write_rows(rows,output);
        std::size_t total_mismatches=0;
        for(const auto& row:rows)
        {
            total_mismatches+=row.mismatches;
        }

        std::cout<<"correctness CSV: "<<output.string()<<'\n';
        std::cout<<"total mismatches: "<<total_mismatches<<'\n';
        return total_mismatches==0?0:2;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"correctness matrix failed: "<<error.what()<<'\n';
        return 1;
    }
}
