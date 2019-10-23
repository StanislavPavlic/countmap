#include <cstdint>
#include <utility>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <future>

#include "index.hpp"
#include "brown_minimizers.hpp"
#include "thread_pool/thread_pool.hpp"

std::vector<minimizer_t> collect_minimizers(const std::unique_ptr<fastaq::FastAQ>& ref, const uint32_t w, const uint32_t k, const uint32_t t) {
  std::vector<minimizer_t> t_minimizers;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool();
  uint32_t start = 0;
  uint32_t end = 0;
  while (end < ref->sequence.size()) {
    while (start < ref->sequence.size() && ref->sequence[end] == 'N') start++;
    if (start == ref->sequence.size()) break;
    end = start;
    while (end < ref->sequence.size() && ref->sequence[end] != 'N') end++;
    uint32_t len = end - start;
    if (len >= w + k + 1) {
      std::vector<std::future<std::vector<minimizer_t>>> thread_futures_ref;
      for (uint32_t tasks = 0; tasks < t - 1; ++tasks) {
        thread_futures_ref.emplace_back(thread_pool->submit(brown::minimizers,
            ref->sequence.c_str() + start + tasks * len / t,
            len / t + w + k - 1,
            k, w));
      }
      thread_futures_ref.emplace_back(thread_pool->submit(brown::minimizers,
            ref->sequence.c_str() + start + (t - 1) * len / t,
            len - (t - 1) * len / t,
            k, w));

      for (uint32_t i = 0; i < t; ++i) {
        thread_futures_ref[i].wait();
        uint32_t offset = start + i * len / t;
        for (auto& el : thread_futures_ref[i].get()) {
          std::get<1>(el) += offset;
          t_minimizers.push_back(el);
        }
      }
    }
    start = end;
  }
  return t_minimizers;
}

// Remove f most frequent minimizers, prepare target minimizers vector for index structure
// Args: t_minimizers - list of target minimizers
//       f            - fraction of most frequent minimizers to be removed
// Return: none
void prep_ref(std::vector<minimizer_t>& t_minimizers, const float f) {
  std::unordered_map<uint64_t, uint32_t> ref_min_frequency;
  for (const auto& minimizer : t_minimizers) {
    ref_min_frequency[std::get<0>(minimizer)]++;
  }

  std::vector<uint32_t> occurences;
  occurences.reserve(ref_min_frequency.size());
  for (const auto& entry : ref_min_frequency) {
    occurences.push_back(entry.second);
  }

  uint32_t position = (uint32_t)((1.0f - f) * (occurences.size() - 1.0f));
  std::sort(occurences.begin(), occurences.end());

  uint32_t cutoff_freq = occurences[position] == 1 ? 2 : occurences[position];

  std::vector<minimizer_t> temp;
  temp.reserve(t_minimizers.size());
  for (const auto& minimizer : t_minimizers) {
    if (ref_min_frequency[std::get<0>(minimizer)] < cutoff_freq) {
      temp.push_back(minimizer);
    }
  }
  std::swap(t_minimizers, temp);
  std::sort(t_minimizers.begin(), t_minimizers.end(),
      [] (const minimizer_t& a, const minimizer_t& b) {
        return (std::get<0>(a) < std::get<0>(b));
      });
  t_minimizers.shrink_to_fit();
}

// Form index structure using minimizers vector
// Args: t_minimizers - list of target minimizers
// Return: map of minimizer values to positions in target minimizers list
std::unordered_map<uint64_t, index_pos_t> index_ref(const std::vector<minimizer_t>& t_minimizers) {
  std::unordered_map<uint64_t, index_pos_t> ref_index;
  uint32_t pos = 0;
  uint32_t num = 0;
  uint64_t prev_min = std::get<0>(t_minimizers[0]);
  for (const auto& minimizer : t_minimizers) {
    if (prev_min != std::get<0>(minimizer)) {
      ref_index[prev_min] = std::make_pair(pos, num);
      pos += num;
      num = 1;
      prev_min = std::get<0>(minimizer); 
    } else {
      num++;
    }
  }
  ref_index[prev_min] = std::make_pair(pos, num);
  return ref_index;
}
