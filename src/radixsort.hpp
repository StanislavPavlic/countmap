#pragma once

// Hit: query position, reference position, relative strand
typedef std::tuple<uint32_t, uint32_t, bool> minimizer_hit_t;

void radixsort(std::vector<minimizer_hit_t>& hits);