#pragma once

#include "mapping_params.hpp"
#include "fastaq.hpp"

// Minimizer: value, position, origin
typedef std::tuple<uint64_t, uint32_t, bool> minimizer_t;
// Index: position, range
typedef std::pair<uint32_t, uint32_t> index_pos_t;
// Hit: query position, reference position, relative strand
typedef std::tuple<uint32_t, uint32_t, bool> minimizer_hit_t;
// Region: two "hits" with swapped positions for reverse complement that "map" the query to the target
typedef std::pair<minimizer_hit_t, minimizer_hit_t> region_t;
// Paired read
typedef std::pair<std::vector<std::unique_ptr<fastaq::FastAQ>>, std::vector<std::unique_ptr<fastaq::FastAQ>>> paired_reads_t;
// Paired read minimizers
typedef std::pair<std::vector<minimizer_t>, std::vector<minimizer_t>> paired_minimizers_t;
// Forward and reverse complement strand hits
typedef std::pair<std::vector<minimizer_hit_t>, std::vector<minimizer_hit_t>> split_hits_t;
// Bin: region, hits
typedef std::unordered_map<uint32_t, std::vector<minimizer_hit_t>> bin_t;
// Paired checked and grouped hits (true candidates)
typedef std::pair<std::vector<std::vector<minimizer_hit_t>>, std::vector<std::vector<minimizer_hit_t>>> paired_checked_t;

void infer_insert_size(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                       const std::vector<minimizer_t>& t_minimizers,
                       const std::unique_ptr<fastaq::FastAQ>& reference, const paired_reads_t& paired_reads,
                       mapping_params_t& parameters);

std::string map_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                       const std::vector<minimizer_t>& t_minimizers,
                       const std::unique_ptr<fastaq::FastAQ>& reference, 
                       const std::vector<std::unique_ptr<fastaq::FastAQ>>& reads, const mapping_params_t& parameters,
                       uint32_t t_start, uint32_t t_end);

std::string map_as_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, const std::unordered_map<uint64_t, index_pos_t>& ref_index_2,
                          const std::vector<minimizer_t>& t_minimizers, const std::vector<minimizer_t>& t_minimizers_2,
                          const std::unique_ptr<fastaq::FastAQ>& reference, const paired_reads_t& paired_reads,
                          const mapping_params_t& parameters, uint32_t t_start, uint32_t t_end);

std::string map_paired(const std::unordered_map<uint64_t, index_pos_t>& ref_index, const std::unordered_map<uint64_t, index_pos_t>& ref_index_2,
                       const std::vector<minimizer_t>& t_minimizers, const std::vector<minimizer_t>& t_minimizers_2,
                       const std::unique_ptr<fastaq::FastAQ>& reference, const paired_reads_t& paired_reads,
                       const mapping_params_t& parameters, uint32_t t_start, uint32_t t_end);