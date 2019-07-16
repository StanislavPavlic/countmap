#include <cstdint>
#include <tuple>
#include <vector>

#include "radixsort.hpp"

// Sort hits by position on the reference
// Args: hits - list of minimizer hits
// Return: none
void radixsort(std::vector<minimizer_hit_t>& hits) {
  if (hits.size() <= 1) {
    return;
  }

  uint32_t max_val = std::get<1>(hits[0]);
  const uint32_t size = hits.size();

  for (uint32_t i = 1; i < size; ++i) {
    if (std::get<1>(hits[i]) > max_val) {
      max_val = std::get<1>(hits[i]);
    }
  }

  std::vector<minimizer_hit_t> out(size);
  for (uint32_t e = 1; max_val / e > 0; e *= 10) {
    uint32_t count[10];

    for (uint32_t i = 0; i < 10; ++i) {
      count[i] = 0;
    }

    for (uint32_t i = 0; i < size; ++i) {
      count[(std::get<1>(hits[i]) / e) % 10]++;
    }

    for (uint32_t i = 1; i < 10; ++i) {
      count[i] += count[i - 1];
    }

    for (uint32_t i = size - 1; i != (uint32_t)(-1); --i) {
      out[count[(std::get<1>(hits[i]) / e) % 10] - 1] = hits[i];
      count[(std::get<1>(hits[i]) / e) % 10]--;
    }
    
    hits.assign(out.begin(), out.end());
  }
}