#pragma once

#include <cstdint>

// Mapping parameters wrapper
struct mapping_params_t {
  bool all;
  int32_t mch;
  int32_t mis;
  int32_t gapo;
  int32_t gape;
  int32_t band;
  uint32_t k;
  uint32_t w;
  uint32_t k_2;
  uint32_t w_2;
  float f;
  uint32_t insert_size;
  float sd;
  uint32_t threshold;
  uint32_t threshold_2;
};