#include <vector>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <algorithm>
#include <iostream>

#include "region.hpp"

// Create hits from matching minimizers on query and target
// Args: fwd          - list of minimizer hits with same relative strand origin
//       rev          - list of minimizer hits with opposite relative strand origin
//       ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       q_minimizers - list of query minimizers
// Return: none
void find_minimizer_hits(std::vector<minimizer_hit_t>& fwd, std::vector<minimizer_hit_t>& rev,
                         const std::unordered_map<uint64_t, index_pos_t>& ref_index,
                         const std::vector<minimizer_t>& t_minimizers, const std::vector<minimizer_t>& q_minimizers) {
  for (const auto& minimizer : q_minimizers) {
    auto found = ref_index.find(std::get<0>(minimizer));
    if (found != ref_index.end()) {
      for (uint32_t i = 0; i < found->second.second; i++) {
        if (std::get<2>(minimizer) == std::get<2>(t_minimizers[found->second.first + i])) {
          fwd.emplace_back(
              std::get<1>(minimizer),
              std::get<1>(t_minimizers[found->second.first + i]),
              0);  
        } else {
          rev.emplace_back(
              std::get<1>(minimizer),
              std::get<1>(t_minimizers[found->second.first + i]),
              1);
        }
      }
    }
  }
}

std::pair<std::vector<minimizer_hit_t>, std::vector<minimizer_hit_t>> circular_check(
    const std::vector<minimizer_hit_t>& hits, const uint32_t region_size, const uint32_t ref_size) {
  std::vector<minimizer_hit_t> circ_s;
  for (uint32_t i = 0; i < hits.size() && std::get<1>(hits[i]) < region_size; ++i) {
    circ_s.push_back(hits[i]);
  }
  std::vector<minimizer_hit_t> circ_e;
  for (uint32_t i = 0; i < hits.size() && std::get<1>(hits[hits.size() - 1 - i]) > ref_size - region_size; ++i) {
    circ_e.push_back(hits[hits.size() - 1 - i]);
  }
  return std::make_pair(circ_s, circ_e);
}

// Count hits by two neighbouring regions, e.g. [0,1], [100, 101], [101, 102], ...
// Args: hits        - list of minimizer hits
//       threshold   - acceptance limit for hit counts in regions
//       region_size - size of region on reference to count hits from
// Return: map of region numbers to hits from that region
bin_t extract_candidates(const std::vector<minimizer_hit_t>& hits, const uint32_t threshold, 
                         const uint32_t region_size) {
  bin_t candidates;
  std::vector<minimizer_hit_t> temp_c;
  std::vector<minimizer_hit_t> temp_n;
  uint32_t current = 0;
  uint32_t next = 0;

  for (uint32_t i = 0; i < hits.size(); ++i) {
    uint32_t reg = std::get<1>(hits[i]) / region_size;
    if (reg == current) {
      temp_c.push_back(hits[i]);
      continue;
    }

    if (reg == current + 1) {
      next = current + 1;
      temp_c.push_back(hits[i]);
      temp_n.push_back(hits[i]);
      continue;
    }

    if (temp_c.size() >= threshold) {
      candidates[current] = std::move(temp_c);
    }

    if (next == current + 1) {
      current = next;
      temp_c = std::move(temp_n);
      temp_n.clear();
      --i;
      continue;
    } else {
      current = reg;
      temp_c.clear();
      temp_c.push_back(hits[i]);
    }
  }

  if (temp_c.size() >= threshold) {
    candidates[current] = std::move(temp_c);
  }
  if (next == current + 1 && temp_n.size() >= threshold) {
    candidates[next] = std::move(temp_n);
  }

  return candidates;
}

// Check paired candidate regions based on read size and insert size and keep viable ones
// E.g. read_size = 100, insert_size = 1000, candidate_region_1 = 10 000 
//        |-> check for candidate_region_2 at 10 000 +/- (900 +/- 1)
// Args: candidates  - pair of maps from region numbers to hits from that region from paired reads
//       insert_size - insert size of fragment from which the reads have been obtained
//       read_size   - read size of second read in pair
//       region_size - size of region on reference from which hits were counted
// Return: pair of corresponding lists containing lists of hits that were checked for insert size
paired_checked_t check_pairing(std::pair<bin_t, bin_t>& candidates, const uint32_t insert_size,
                               const uint32_t read_size, const uint32_t ref_size) {
  paired_checked_t checked;
  bool circd = false;
  minimizer_hit_t flag = {-1, -1, 0};

  if (candidates.first.find(0) != candidates.first.end() && candidates.first.find(ref_size / read_size) != candidates.first.end()) {
    bin_t::const_iterator found;
    // std::cerr << "ENTER 1" << std::endl;

    if (!std::get<2>(candidates.first[0][0])) {
      if ((found = candidates.second.find((insert_size - read_size) / read_size)) != candidates.second.end()
          || (found = candidates.second.find(((insert_size - read_size) / read_size) - 1)) != candidates.second.end()
          || (found = candidates.second.find(((insert_size - read_size) / read_size) + 1)) != candidates.second.end()) {
        // std::cerr << std::get<1>(candidates.first[0][0]) << ", " << std::get<1>(candidates.first[ref_size / read_size][0]) << ", " << std::get<1>(found->second[0]) << std::endl;
        checked.first.emplace_back(candidates.first[0].begin(), candidates.first[0].end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
        checked.first.emplace_back(candidates.first[ref_size / read_size].begin(), candidates.first[ref_size / read_size].end());
        checked.second.emplace_back(1, flag);
        circd = true;
      }
    } else {
      if ((found = candidates.second.find((ref_size - (insert_size - read_size)) / read_size)) != candidates.second.end()
          || (found = candidates.second.find((ref_size - ((insert_size - read_size)) / read_size) - 1)) != candidates.second.end()
          || (found = candidates.second.find((ref_size - ((insert_size - read_size)) / read_size) + 1)) != candidates.second.end()) {
        // std::cerr << std::get<1>(candidates.first[0][0]) << ", " << std::get<1>(candidates.first[ref_size / read_size][0]) << ", " << std::get<1>(found->second[0]) << std::endl;
        checked.first.emplace_back(candidates.first[0].begin(), candidates.first[0].end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
        checked.first.emplace_back(candidates.first[ref_size / read_size].begin(), candidates.first[ref_size / read_size].end());
        checked.second.emplace_back(1, flag);
        circd = true;
      }
    }
    // std::cerr << "EXIT 1" << std::endl;
  }
  if (candidates.second.find(0) != candidates.second.end() && candidates.second.find(ref_size / read_size) != candidates.second.end()) {
    bin_t::const_iterator found;
    // std::cerr << "ENTER 2" << std::endl;

    if (!std::get<2>(candidates.second[0][0])) {
      if ((found = candidates.first.find((insert_size - read_size) / read_size)) != candidates.first.end()
          || (found = candidates.first.find(((insert_size - read_size) / read_size) - 1)) != candidates.first.end()
          || (found = candidates.first.find(((insert_size - read_size) / read_size) + 1)) != candidates.first.end()) {
        // std::cerr << std::get<1>(candidates.second[0][0]) << ", " << std::get<1>(candidates.second[ref_size / read_size][0]) << ", " << std::get<1>(found->second[0]) << std::endl;
        checked.first.emplace_back(found->second.begin(), found->second.end());
        checked.second.emplace_back(candidates.second[0].begin(), candidates.second[0].end());
        checked.first.emplace_back(1, flag);
        checked.second.emplace_back(candidates.second[ref_size / read_size].begin(), candidates.second[ref_size / read_size].end());
        circd = true;
      }
    } else {
      if ((found = candidates.first.find((ref_size - (insert_size - read_size)) / read_size)) != candidates.first.end()
          || (found = candidates.first.find((ref_size - ((insert_size - read_size)) / read_size) - 1)) != candidates.first.end()
          || (found = candidates.first.find((ref_size - ((insert_size - read_size)) / read_size) + 1)) != candidates.first.end()) {
        // std::cerr << std::get<1>(candidates.second[0][0]) << ", " << std::get<1>(candidates.second[ref_size / read_size][0]) << ", " << std::get<1>(found->second[0]) << std::endl;
        checked.first.emplace_back(found->second.begin(), found->second.end());
        checked.second.emplace_back(candidates.second[0].begin(), candidates.second[0].end());
        checked.first.emplace_back(1, flag);
        checked.second.emplace_back(candidates.second[ref_size / read_size].begin(), candidates.second[ref_size / read_size].end());
        circd = true;
      }
    }
    // std::cerr << "EXIT 2" << std::endl;
  }
  if (circd) {
    return checked;
  }

  for (const auto& bin : candidates.first) {
    bin_t::const_iterator found;

    if (!std::get<2>(bin.second[0])) {
      if ((found = candidates.second.find(bin.first + ((insert_size - read_size) / read_size))) != candidates.second.end()
          || (found = candidates.second.find(bin.first + ((insert_size - read_size) / read_size) - 1)) != candidates.second.end()
          || (found = candidates.second.find(bin.first + ((insert_size - read_size) / read_size) + 1)) != candidates.second.end()) {
        checked.first.emplace_back(bin.second.begin(), bin.second.end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
      }
    } else {
      if ((found = candidates.second.find(bin.first - ((insert_size - read_size) / read_size))) != candidates.second.end()
          || (found = candidates.second.find(bin.first - ((insert_size - read_size) / read_size) - 1)) != candidates.second.end()
          || (found = candidates.second.find(bin.first - ((insert_size - read_size) / read_size) + 1)) != candidates.second.end()) {
        checked.first.emplace_back(bin.second.begin(), bin.second.end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
      }
    }
  }

  return checked;
}

// Perform LIS algorithm on reference positions and form region
// Args: hits - list of minimizer hits
// Return: region on query and target
region_t find_region(std::vector<minimizer_hit_t>& hits) {
  if (hits.empty()) {
    return std::make_pair(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0));
  }
  
  std::sort(hits.begin(), hits.end(), 
            [] (const minimizer_hit_t& a, const minimizer_hit_t& b) {
              if (std::get<2>(a) == 0) {
                if (std::get<0>(a) == std::get<0>(b)) {
                  return std::get<1>(a) < std::get<1>(b);
                }
                return std::get<0>(a) < std::get<0>(b);
              } else {
                if (std::get<0>(a) == std::get<0>(b)) {
                  return std::get<1>(a) < std::get<1>(b);
                }
                return std::get<0>(a) > std::get<0>(b);
              }
            }
  );

  std::vector<region_t> lis(hits.size(),
      std::make_pair(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0)));

  uint32_t len = 1;
  lis[0] = std::make_pair(hits[0], hits[0]);

  for (uint32_t i = 1; i < hits.size(); ++i) {
    if (std::get<1>(hits[i]) > std::get<1>(lis[len - 1].second)) {
      lis[len] = lis[len - 1];
      lis[len].second = hits[i];
      len++;
    } else {
      auto pair_it = std::upper_bound(lis.begin(), lis.begin() + len, hits[i],
          [] (const minimizer_hit_t& a,
              const region_t& b) {
            return (std::get<1>(a) < std::get<1>(b.second));
          });
      if (pair_it == lis.begin()) {
        *pair_it = std::make_pair(hits[i], hits[i]);
      } else {
        *pair_it = std::make_pair((pair_it - 1)->first, hits[i]);
      }
    }
  }

  region_t region = lis[len - 1];

  if (std::get<2>(region.first) == 1) {
    uint32_t temp = std::get<0>(region.first);
    std::get<0>(region.first) = std::get<0>(region.second);
    std::get<0>(region.second) = temp;
  }

  return region; 
}

// Expand region to the whole read size
// Args: reg       - region to be expanded
//       read_size - size of read to which the region is expanded
//       k         - length of k-mers used to obtain minimizers
//       max_size  - size of target
void expand_region(region_t& reg, uint32_t read_size, uint32_t k, uint32_t max_size) {
  if (std::get<2>(reg.first) == 0) {
    if (std::get<1>(reg.first) < std::get<0>(reg.first)) {
      std::get<0>(reg.first) -= std::get<1>(reg.first);
      std::get<1>(reg.first) = 0;
    } else {
      std::get<1>(reg.first) -= std::get<0>(reg.first);
      std::get<0>(reg.first) = 0;
    }
    if (max_size < std::get<1>(reg.second) + read_size - std::get<0>(reg.second)) {
      std::get<0>(reg.second) += max_size - std::get<1>(reg.second);
      std::get<1>(reg.second) = max_size;
    } else {
      std::get<1>(reg.second) += read_size - std::get<0>(reg.second);
      std::get<0>(reg.second) = read_size;
    }
  } else {
    if (std::get<1>(reg.first) < read_size - (std::get<0>(reg.second) + k)) {
      std::get<0>(reg.second) += std::get<1>(reg.first) + k;
      std::get<1>(reg.first) = 0;
    } else {
      std::get<1>(reg.first) -= read_size - (std::get<0>(reg.second) + k);
      std::get<0>(reg.second) = read_size;
    }
    if (max_size < std::get<1>(reg.second) + (std::get<0>(reg.first) + k)) {
      std::get<0>(reg.first) -= max_size - std::get<1>(reg.second) - k;
      std::get<1>(reg.second) = max_size;
    } else {
      std::get<1>(reg.second) += std::get<0>(reg.first) + k;
      std::get<0>(reg.first) = 0;
    }
  }
}