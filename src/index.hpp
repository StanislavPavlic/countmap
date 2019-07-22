#pragma once

#include "fastaq.hpp"
#include "mapping_params.hpp"

// Minimizer: value, position, origin
typedef std::tuple<uint64_t, uint32_t, bool> minimizer_t;
// Index: position, range
typedef std::pair<uint32_t, uint32_t> index_pos_t;

std::vector<minimizer_t> collect_minimizers(const std::unique_ptr<fastaq::FastAQ>& ref, const mapping_params_t& parameters, const uint32_t t);

void prep_ref(std::vector<minimizer_t>& t_minimizers, const float f);

std::unordered_map<uint64_t, index_pos_t> index_ref(const std::vector<minimizer_t>& t_minimizers);