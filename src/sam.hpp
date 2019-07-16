#pragma once

#include "mapping_params.hpp"
#include "mapping.hpp"

// Hit: query position, reference position, relative strand
typedef std::tuple<uint32_t, uint32_t, bool> minimizer_hit_t;
// Region: two "hits" with swapped positions for reverse complement that "map" the query to the target
typedef std::pair<minimizer_hit_t, minimizer_hit_t> region_t;

std::string sam_format(const mapping_t& m);

mapping_t single_mapping(const std::string& qname, const std::string& query, 
                         const std::string& qual, const std::string& rname, 
                         const std::string& ref, const region_t& region, 
                         const mapping_params_t& parameters, 
                         const std::pair<int32_t, int32_t>& clipped);

std::pair<mapping_t, mapping_t> pair_mapping(const std::string& qname,
    const std::string& query1, const std::string& qual1, const std::string& query2, const std::string& qual2,
    const std::string& rname, const std::string& ref,
    const std::pair<region_t, region_t>& region_pair, const mapping_params_t& parameters,
    const std::pair<int32_t, int32_t>& clipped1, const std::pair<int32_t, int32_t>& clipped2);

std::string unmapped_sam(const std::string& qname, const std::string& query, const std::string& qual, 
                         bool pair, bool first, bool last);