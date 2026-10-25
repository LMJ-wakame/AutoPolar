#include "polar_scl/construction.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>

namespace polar_scl
{
    namespace
    {
        bool is_power_of_two(std::size_t x)
        {
            return x != 0 && (x & (x - 1)) == 0;
        }

        std::vector<Bit> polynomial_remainder(
        const std::vector<Bit>& bits,
        std::size_t width,
        std::uint32_t polynomial)
        {
            if(width < 1 || width > 31) throw std::invalid_argument("CRC width must lie in [1,31]");
            if(polynomial >= (std::uint32_t{1} << width))
            throw std::invalid_argument("polynomial must omit the leading term and fit CRC width");
            for(Bit b : bits) if(b > 1) throw std::invalid_argument("CRC input must be binary");

            std::vector<Bit> divisor(width + 1, 0);
            divisor[0] = 1;
            for(std::size_t i = 0; i < width; ++i)
            {
                const std::size_t shift = width - 1 - i;
                divisor[i + 1] = static_cast<Bit>((polynomial >> shift) & 1U);
            }

            std::vector<Bit> work = bits;
            if(work.size() < width) work.insert(work.begin(), width - work.size(), Bit{0}

            );
            for(std::size_t pos = 0; pos + width < work.size(); ++pos)
            {
                if(work[pos] == 0) continue;
                for(std::size_t j = 0; j <= width; ++j) work[pos + j] ^= divisor[j];
            }

            if(work.size() < width) return std::vector<Bit>(width, 0);
            return std::vector<Bit>(work.end() - static_cast<std::ptrdiff_t>(width), work.end());
        }
    }

    // namespace

    std::vector<std::size_t> polarization_weight_order(std::size_t length, double beta)
    {
        if(!is_power_of_two(length)) throw std::invalid_argument("length must be a positive power of two");
        if(beta <= 1.0) throw std::invalid_argument("beta must be greater than one");
        std::size_t stages = 0;
        for(auto n = length; n > 1; n >>= 1) ++stages;
        std::vector<std::size_t> order(length);
        std::iota(order.begin(), order.end(), 0);
        auto weight = [&](std::size_t index)
        {
            double result = 0.0;
            for(std::size_t stage = 0; stage < stages; ++stage)
            if(((index >> stage) & 1ULL) != 0ULL) result += std::pow(beta, static_cast<double>(stage));
            return result;
        };
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
        {
            const auto wa = weight(a), wb = weight(b);
            return wa < wb || (wa == wb && a < b);
        }

        );
        return order;
    }

    std::vector<bool> pw_frozen_pattern(std::size_t length, std::size_t information_bits)
    {
        if(information_bits > length) throw std::invalid_argument("information_bits must lie in [0,length]");
        auto order = polarization_weight_order(length);
        std::set<std::size_t> info(order.end() - static_cast<std::ptrdiff_t>(information_bits), order.end());
        std::vector<bool> frozen(length, true);
        for(std::size_t i = 0; i < length; ++i) frozen[i] = !info.contains(i);
        return frozen;
    }

    std::vector<Bit> place_information(const std::vector<Bit>& information, const std::vector<bool>& frozen)
    {
        for(Bit bit : information) if(bit > 1) throw std::invalid_argument("information must be binary");
        const auto count = static_cast<std::size_t>(std::count(frozen.begin(), frozen.end(), false));
        if(information.size() != count) throw std::invalid_argument("information length does not match non-frozen positions");
        std::vector<Bit> u(frozen.size(), 0);
        std::size_t cursor = 0;
        for(std::size_t i = 0; i < frozen.size(); ++i) if(!frozen[i]) u[i] = information[cursor++];
        return u;
    }

    std::vector<Bit> extract_information(const std::vector<Bit>& bits, const std::vector<bool>& frozen)
    {
        if(bits.size() != frozen.size()) throw std::invalid_argument("bits and frozen pattern must have equal length");
        std::vector<Bit> result;
        for(std::size_t i = 0; i < bits.size(); ++i) if(!frozen[i]) result.push_back(bits[i]);
        return result;
    }

    std::vector<Bit> append_crc(const std::vector<Bit>& payload, std::size_t width, std::uint32_t polynomial)
    {
        std::vector<Bit> padded = payload;
        padded.insert(padded.end(), width, 0);
        auto remainder = polynomial_remainder(padded, width, polynomial);
        auto result = payload;
        result.insert(result.end(), remainder.begin(), remainder.end());
        return result;
    }

    bool crc_ok(const std::vector<Bit>& information, std::size_t width, std::uint32_t polynomial)
    {
        if(information.size() < width) return false;
        const auto remainder = polynomial_remainder(information, width, polynomial);
        return std::all_of(remainder.begin(), remainder.end(), [](Bit b)
        {
            return b == 0;
        }

        );
    }

    std::pair<std::vector<Bit>, std::vector<Bit>> encode_payload(
    const std::vector<Bit>& payload,
    const std::vector<bool>& frozen,
    std::size_t crc_width,
    std::uint32_t polynomial)
    {
        auto info = append_crc(payload, crc_width, polynomial);
        auto u = place_information(info, frozen);
        return {u, polar_transform(u)};
    }

    DecoderPath select_crc_path(
    const std::vector<DecoderPath>& paths,
    const std::vector<bool>& frozen,
    std::size_t crc_width,
    std::uint32_t polynomial)
    {
        if(paths.empty()) throw std::invalid_argument("at least one decoder path is required");
        auto sorted = paths;
        std::sort(sorted.begin(), sorted.end(), decoder_path_less);
        for(const auto& path : sorted) if(crc_ok(extract_information(path.u, frozen), crc_width, polynomial)) return path;
        return sorted.front();
    }
}

// namespace polar_scl
