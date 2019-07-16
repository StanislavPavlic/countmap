#pragma once

// Minimizer: value, position, origin
typedef std::tuple<uint64_t, uint32_t, bool> minimizer_t;
// Index: position, range
typedef std::pair<uint32_t, uint32_t> index_pos_t;
// Hit: query position, reference position, relative strand
typedef std::tuple<uint32_t, uint32_t, bool> minimizer_hit_t;
// Region: two "hits" with swapped positions for reverse complement that "map" the query to the target
typedef std::pair<minimizer_hit_t, minimizer_hit_t> region_t;
// Bin: region, hits
typedef std::unordered_map<uint32_t, std::vector<minimizer_hit_t>> bin_t;
// Paired checked and grouped hits (true candidates)
typedef std::pair<std::vector<std::vector<minimizer_hit_t>>, std::vector<std::vector<minimizer_hit_t>>> paired_checked_t;

void find_minimizer_hits(std::vector<minimizer_hit_t>& fwd, std::vector<minimizer_hit_t>& rev,
                         const std::unordered_map<uint64_t, index_pos_t>& ref_index,
                         const std::vector<minimizer_t>& t_minimizers, const std::vector<minimizer_t>& q_minimizers);

bin_t extract_candidates(const std::vector<minimizer_hit_t>& hits, const uint32_t threshold, 
                         const uint32_t region_size);

paired_checked_t check_pairing(std::pair<bin_t, bin_t>& candidates, const uint32_t insert_size,
                               const uint32_t read_size);

region_t find_region(std::vector<minimizer_hit_t>& hits);

void expand_region(region_t& reg, uint32_t read_size, uint32_t k, uint32_t max_size);