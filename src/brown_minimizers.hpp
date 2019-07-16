#pragma once

namespace brown {

std::vector<std::tuple<uint64_t, uint32_t, bool>> minimizers(
    const char* sequence,
    uint32_t sequence_length,
    uint32_t k,
    uint32_t window_length);

}