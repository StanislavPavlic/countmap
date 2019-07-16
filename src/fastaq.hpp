#pragma once

#include <string>
#include <vector>
#include <limits>
#include <cstdint>

#include "bioparser/bioparser.hpp"

namespace fastaq {

typedef struct {
  uint32_t num;
  double avg;
  uint32_t max;
  uint32_t min;
} stats;
  
class FastAQ {
public:
  std::string name;
  std::string sequence;
  std::string quality;

  FastAQ(
    const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length) : FastAQ(name, name_length, sequence, sequence_length, "*", 0){}

  FastAQ(
    const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length) {
      this->name = {name, name_length};
      this->sequence = {sequence, sequence_length};
      this->quality = {quality, quality_length};
  }

  static stats print_statistics(
    const std::vector<std::unique_ptr<FastAQ>> &fastaq_objects,
    const std::string file) {
      uint32_t num = fastaq_objects.size();
      double average = 0;
      uint32_t max = 0;
      uint32_t min = std::numeric_limits<uint32_t>::max();
      for (uint32_t i = 0; i < num; i++) {
        average += fastaq_objects[i]->sequence.size();
        if (fastaq_objects[i]->sequence.size() > max) {
          max = fastaq_objects[i]->sequence.size();
        }
        if (fastaq_objects[i]->sequence.size() < min) {
          min = fastaq_objects[i]->sequence.size();
        }
      }
      average /= num;
      fprintf(stderr, "[countmap-stats] stats for: %s\n"
                      "  number of sequences: %u\n"
                      "  average length:      %g\n"
                      "  maximum length:      %u\n"
                      "  minimum length:      %u\n",
                      file.c_str(), num, average, max, min);
      return {num, average, max, min};
  }
};

}