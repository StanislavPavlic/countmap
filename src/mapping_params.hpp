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
  float f;
  uint32_t insert_size;
  float sd;
  uint32_t threshold;
};