#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <unordered_set>
#include <tuple>
#include <queue>
#include <limits>
#include <algorithm>
#include <functional>
#include <map>

#include "brown_minimizers.hpp"

namespace brown {

std::map<char, uint64_t> char_to_val = {{'C', 0}, {'A', 1}, {'T', 2}, {'U', 2}, {'G', 3}};

typedef std::tuple<uint64_t, uint32_t, bool> triplet_t;

struct triplet_hash {
  std::size_t operator() (const triplet_t& k) const noexcept {
    return std::hash<uint64_t>{}(std::get<0>(k)) ^ std::hash<uint32_t>{}(std::get<1>(k)) ^ std::hash<bool>{}(std::get<2>(k));
  }
};

struct triplet_equal {
  bool operator() (const triplet_t& v0, const triplet_t& v1) const noexcept {
    return (
            std::get<0>(v0) == std::get<0>(v1) &&
            std::get<1>(v0) == std::get<1>(v1) &&
            std::get<2>(v0) == std::get<2>(v1)
           );
  }
};

struct triplet_ordering {
  bool operator() (const triplet_t& v0, const triplet_t& v1) const noexcept {
    return (std::get<1>(v0) < std::get<1>(v1));
  }
};

inline uint64_t value(const char* sequence, uint32_t pos, uint32_t k) {
  uint64_t val = 0;
  for (uint32_t i = 0; i < k; ++i) {
    val = (val << 2) | char_to_val[sequence[pos + i]];
  }
  return val;
}

inline uint64_t value_reverse_complement(const char* sequence,
                                         uint32_t pos,
                                         uint32_t k) {

  uint64_t val = 0;
  for (uint32_t i = k - 1; i < (uint32_t)(-1); --i) {
    val = (val << 2) | (~char_to_val[sequence[pos + i]] & 3);
  }
  return val;
}

void interior_minimizers_fill(
    std::unordered_set<triplet_t, triplet_hash, triplet_equal>& minimizers_set,
    const char* sequence,
    uint32_t sequence_length,
    uint32_t k,
    uint32_t window_length) {

  uint64_t mask = ((uint64_t)1 << (2 * k)) - 1;
  uint64_t original;
  uint64_t rev_com;
  uint64_t foriginal;
  uint64_t frev_com;

  bool cache_empty = true;
  triplet_t cache;
  for (uint32_t i = 0; i <= sequence_length - window_length - k + 1; ++i) {
    uint64_t m = std::numeric_limits<uint64_t>::max();

    if (cache_empty) {
      foriginal = value(sequence, i, k) >> 2;
      frev_com = (value_reverse_complement(sequence, i, k) << 2) & mask;
    } else {
      foriginal = original;
      frev_com = rev_com;
    }
    original = foriginal;
    rev_com = frev_com;

    uint32_t start;
    start = cache_empty ? 0 : window_length - 1;

    for (uint32_t j = start; j < window_length; ++j) {
      uint32_t pos = k - 1 + i + j;
      original = ((original << 2) | char_to_val[sequence[pos]]) & mask;
      rev_com = ((rev_com >> 2) | ((~char_to_val[sequence[pos]] & 3) << (2 * (k - 1)))) & mask;
      if (original != rev_com) {
        m = std::min(m, std::min(original, rev_com));
      }
    }

    if (!cache_empty) {
      if (std::get<1>(cache) >= i && std::get<0>(cache) <= m) {
        m = std::get<0>(cache);
      } else {
        cache_empty = true;
      }
    }

    original = foriginal;
    rev_com = frev_com;

    for (uint32_t j = start; j < window_length; ++j) {
      uint32_t pos = k - 1 + i + j;
      original = ((original << 2) | char_to_val[sequence[pos]]) & mask;
      rev_com = ((rev_com >> 2) | ((~char_to_val[sequence[pos]] & 3) << (2 * (k - 1)))) & mask;
      if (original < rev_com && original == m) {
        minimizers_set.emplace(m, i + j, 0);
        cache = std::make_tuple(m, i + j, 0);
        cache_empty = false;
      } else if (rev_com < original && rev_com == m) {
        minimizers_set.emplace(m, i + j, 1);
        cache = std::make_tuple(m, i + j, 1);
        cache_empty = false;
      }
    }
  }
}

void end_minimizers_fill(
    std::unordered_set<triplet_t, triplet_hash, triplet_equal>& minimizers_set,
    uint32_t k,
    uint32_t window_length,
    std::function<uint32_t(uint32_t)> position,
    std::function<uint64_t(uint32_t, uint64_t)> start_o,
    std::function<uint64_t(uint32_t, uint64_t)> start_r,
    std::function<uint64_t(uint64_t, int, uint64_t)> insert_o,
    std::function<uint64_t(uint64_t, int, uint64_t)> insert_r) {

  uint64_t mask = ((uint64_t)1 << (2 * k)) - 1;
  uint64_t original;
  uint64_t rev_com;
  uint64_t m = std::numeric_limits<uint64_t>::max();

  for (uint32_t i = 0; i < window_length - 1; ++i) {
    if (i == 0) {
      original = start_o(position(i), mask);
      rev_com = start_r(position(i), mask);
    }

    original = insert_o(original, position(i), mask);
    rev_com = insert_r(rev_com, position(i), mask);

    if (original != rev_com) {
      m = std::min(m, std::min(original, rev_com));
    }

    if (original < rev_com && original == m) {
      minimizers_set.emplace(m, position(i), 0);
    } else if (rev_com < original && rev_com == m) {
      minimizers_set.emplace(m, position(i), 1);
    }
  }
}

std::vector<triplet_t> minimizers(
    const char* sequence,
    uint32_t sequence_length,
    uint32_t k,
    uint32_t window_length) {

  if (k > 32) {
    fprintf(stderr, "[brown::minimizers] error: Largest supported value for k is 32.\n"
                    "  k = %d.\n", k);
    exit(1);
  }
  if (sequence_length < window_length + k - 1) {
    fprintf(stderr, "[brown::minimizers] error: Sequence length too short for given parameters.\n"
                    "  Length = %d, k = %d, window length = %d.\n",
        sequence_length, k, window_length);
    exit(1);
  }

  char_to_val['N'] = std::rand() % 4;

  std::unordered_set<triplet_t, triplet_hash, triplet_equal> minimizers_set;
  
  interior_minimizers_fill(minimizers_set, sequence, sequence_length, k, window_length);

  end_minimizers_fill(minimizers_set, k, window_length,
                      [] (uint32_t i) {return i;},
                      [&k, &sequence] (uint32_t pos, uint64_t mask) {
                        return (value(sequence, pos, k) >> 2) & mask;
                      },
                      [&k, &sequence] (uint32_t pos, uint64_t mask) {
                        return (value_reverse_complement(sequence, pos, k) << 2) & mask;
                      },
                      [&k, &sequence] (uint64_t original, uint32_t pos, uint64_t mask) {
                        return ((original << 2) | char_to_val[sequence[k - 1 + pos]]) & mask;
                      },
                      [&k, &sequence] (uint64_t rev_com, uint32_t pos, uint64_t mask) {
                        return ((rev_com >> 2) | ((~char_to_val[sequence[k - 1 + pos]] & 3) << (2 * (k - 1)))) & mask;
                      }
                     );

  end_minimizers_fill(minimizers_set, k, window_length,
                      [&sequence_length, &k] (uint32_t i) {return sequence_length - k - i;},
                      [&k, &sequence] (uint32_t pos, uint64_t mask) {
                        return (value(sequence, pos, k) << 2) & mask;
                      },
                      [&k, &sequence] (uint32_t pos, uint64_t mask) {
                        return (value_reverse_complement(sequence, pos, k) >> 2) & mask;
                      },
                      [&k, &sequence] (uint64_t original, uint32_t pos, uint64_t mask) {
                        return ((original >> 2) | (char_to_val[sequence[pos]] << (2 * (k - 1)))) & mask;
                      },
                      [&sequence] (uint64_t rev_com, uint32_t pos, uint64_t mask) {
                        return ((rev_com << 2) | (~char_to_val[sequence[pos]] & 3)) & mask;
                      }
                     );

  std::vector<triplet_t> minimizers_vector(
      minimizers_set.begin(), minimizers_set.end());

  return minimizers_vector;
}

}