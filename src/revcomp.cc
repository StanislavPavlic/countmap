#include <cstdint>
#include <string>
#include <unordered_map>

#include "revcomp.hpp"

std::unordered_map<char, char> complement_map = {{'C', 'G'}, {'A', 'T'}, {'T', 'A'}, {'U', 'A'}, {'G', 'C'}, {'N', 'N'}};

// Create reverse complement of sequence
// Args: original - sequence whose subsequence will be reversed and complemented
//       pos      - starting position of subsequence on original
//       len      - length of subsequence to be reversed and complemented
// Return: reverse complement of original sequence subsequence
std::string reverse_complement(const std::string& original, uint32_t pos, uint32_t length) {
  std::string rc(original.begin() + pos, original.begin() + pos + length);
  uint32_t j = pos + length - 1;

  for (uint32_t i = 0; i < length; ++i) {
    rc[i] = complement_map[original[j--]];
  }
  
  return rc;
}