#pragma once

#include "polar_scl/decoder.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace polar_scl
{
    std::vector<std::size_t> polarization_weight_order(
    std::size_t length,
    double beta = 1.1892071150027210667);

    std::vector<bool> pw_frozen_pattern(
    std::size_t length,
    std::size_t information_bits);

    std::vector<Bit> place_information(
    const std::vector<Bit>& information,
    const std::vector<bool>& frozen);

    std::vector<Bit> extract_information(
    const std::vector<Bit>& bits,
    const std::vector<bool>& frozen);

    std::vector<Bit> append_crc(
    const std::vector<Bit>& payload,
    std::size_t width = 8,
    std::uint32_t polynomial = 0x07);

    bool crc_ok(
    const std::vector<Bit>& information,
    std::size_t width = 8,
    std::uint32_t polynomial = 0x07);

    std::pair<std::vector<Bit>, std::vector<Bit>> encode_payload(
    const std::vector<Bit>& payload,
    const std::vector<bool>& frozen,
    std::size_t crc_width = 8,
    std::uint32_t polynomial = 0x07);

    DecoderPath select_crc_path(
    const std::vector<DecoderPath>& paths,
    const std::vector<bool>& frozen,
    std::size_t crc_width = 8,
    std::uint32_t polynomial = 0x07);
}

// namespace polar_scl
