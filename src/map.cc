#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include "map.hpp"
#include "mapping.hpp"
#include "radixsort.hpp"
#include "revcomp.hpp"
#include "region.hpp"
#include "brown_minimizers.hpp"
#include "sam.hpp"

// Clip ends of sequence that contain multiple Ns
// Args: seq - sequence to be clipped
//       k   - k-mer length
//       w   - window length
// Return: size of clippings; if negative, sequence should be rejected
std::pair<int32_t, int32_t> clip(const std::string& seq, const uint32_t k, const uint32_t w) {
  std::size_t found1 = seq.find_first_of("ACGT");
  if (found1 == std::string::npos || seq.size() - found1 < w + k - 1) {
    return std::make_pair(-1, -1);
  }
  std::size_t found2 = seq.find("NNN", found1);
  if (found2 == std::string::npos) {
    return std::make_pair(found1, 0);
  }
  if (found2 - found1 < w + k - 1) {
    return std::make_pair(-1, -1);
  }
  return std::make_pair(found1, seq.size() - found2);
}

// Process (as) single-end sequencing read
// Args: sam          - SAM format output string
//       ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       read         - FastAQ representation of read
//       parameters   - mapping parameters
// Return: mappings
std::vector<mapping_t> process_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                                      const std::vector<minimizer_t>& t_minimizers, const std::unique_ptr<fastaq::FastAQ>& reference, 
                                      const std::unique_ptr<fastaq::FastAQ>& read, const mapping_params_t& parameters) {
  std::vector<mapping_t> mappings;

  std::pair<int32_t, int32_t> clipped = clip(read->sequence, parameters.k, parameters.w);
  if (clipped.first < 0) {
    return mappings;
  }

  std::vector<minimizer_t> q_minimizers = brown::minimizers(read->sequence.c_str() + clipped.first, 
                                                            read->sequence.size() - clipped.first - clipped.second,
                                                            parameters.k, parameters.w);
  if (clipped.first) {
    for (auto& min : q_minimizers) {
      std::get<1>(min) += clipped.first;
    }
  }

  split_hits_t hits;
  find_minimizer_hits(hits.first, hits.second, ref_index, t_minimizers, q_minimizers);

  radixsort(hits.first);
  radixsort(hits.second);

  std::pair<bin_t, bin_t> candidates(extract_candidates(hits.first, parameters.threshold, read->sequence.size()),
                                     extract_candidates(hits.second, parameters.threshold, read->sequence.size()));

  if (candidates.first.size() == 0 && candidates.second.size() == 0) {
    return mappings;
  }

  std::unordered_set<uint32_t> processed;
  for (auto& bin : candidates.first) {
    region_t reg = find_region(bin.second);
    expand_region(reg, read->sequence.size(), parameters.k, reference->sequence.size());

    if (processed.find(std::get<1>(reg.first)) != processed.end()) {
      continue;
    }
    processed.insert(std::get<1>(reg.first));

    mappings.emplace_back(single_mapping(read->name, read->sequence, read->quality,
                                         reference->name, reference->sequence, reg, 
                                         parameters, clipped));
  }
  for (auto& bin : candidates.second) {
    region_t reg = find_region(bin.second);
    expand_region(reg, read->sequence.size(), parameters.k, reference->sequence.size());

    if (processed.find(std::get<1>(reg.first)) != processed.end()) {
      continue;
    }
    processed.insert(std::get<1>(reg.first));

    std::string rc = reverse_complement(read->sequence, 0, read->sequence.size());
    std::string rq = std::string(read->quality.rbegin(), read->quality.rend());

    mappings.emplace_back(single_mapping(read->name, rc, rq,
                                         reference->name, reference->sequence, reg, 
                                         parameters, clipped));
  }

  std::sort(mappings.begin(), mappings.end(), 
            [] (const mapping_t& a, const mapping_t& b) {
              return a.as > b.as;
            }
  );

  if (parameters.all) {
    for (uint32_t i = 0; i < mappings.size(); ++i) {
      mappings[i].mapq /= mappings.size();
      if (i) {
        mappings[i].flag |= 0x100;
        mappings[i].mapq = 0;
      }
    }
  } else {
    uint32_t i = 1;
    int32_t as = mappings[0].as;
    while (i < mappings.size() && mappings[i].as == as) ++i;
    uint32_t pick = std::rand() % i;
    mappings[0] = mappings[pick];
    mappings[0].mapq /= mappings.size();
    mappings.erase(mappings.begin() + 1, mappings.end());
  }

  return mappings;
}

// Process paired-end sequencing reads
// Args: mappings   - list of SAM format output string with priority values
//       checked    - pair of corresponding lists containing lists of hits that were checked for insert size
//       reference  - FastAQ representation of reference
//       first      - FastAQ representation of first read from pair
//       second     - FastAQ representation of second read from pair
//       parameters - mapping parameters
// Return: mappings
void process_pairs(std::vector<std::pair<mapping_t, mapping_t>>& mappings,
    paired_checked_t& checked, const std::unique_ptr<fastaq::FastAQ>& reference, 
    const std::unique_ptr<fastaq::FastAQ>& first, const std::unique_ptr<fastaq::FastAQ>& second,
    const mapping_params_t& parameters, 
    const std::pair<int32_t, int32_t>& clipped1, 
    const std::pair<int32_t, int32_t>& clipped2) {
  std::unordered_set<uint32_t> processed;

  for (uint32_t j = 0; j < checked.first.size(); ++j) {
    std::pair<region_t, region_t> region_pair(find_region(checked.first[j]),
                                              find_region(checked.second[j]));

    expand_region(region_pair.first, first->sequence.size(), 
                  parameters.k, reference->sequence.size());
    expand_region(region_pair.second, second->sequence.size(), 
                  parameters.k, reference->sequence.size());

    if (processed.find(std::get<1>(region_pair.first.first)) != processed.end()) continue;
    processed.insert(std::get<1>(region_pair.first.first));

    std::string rc;
    std::string rq;

    std::string* query1 = std::get<2>(region_pair.first.first)
                          ? &(rc = reverse_complement(first->sequence, 0, first->sequence.size()))
                          : &(first->sequence);
    std::string* quality1 = std::get<2>(region_pair.first.first)
                            ? &(rq = std::string(first->quality.rbegin(), first->quality.rend()))
                            : &(first->quality);
    std::string* query2 = std::get<2>(region_pair.second.first)
                          ? &(rc = reverse_complement(second->sequence, 0, second->sequence.size()))
                          : &(second->sequence);
    std::string* quality2 = std::get<2>(region_pair.second.first)
                            ? &(rq = std::string(second->quality.rbegin(), second->quality.rend()))
                            : &(second->quality);

    mappings.emplace_back(
      pair_mapping(
        first->name,
        *query1, *quality1, 
        *query2, *quality2, 
        reference->name, reference->sequence, 
        region_pair, parameters, clipped1, clipped2)
    );
  }
}

// Infer insert size based on relatively good mappings of paired-end reads mapped
// as single-end
// Args: ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       paired_reads - pair of lists of FastAQ representations of paired reads
//       parameters   - mapping parameters
void infer_insert_size(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                       const std::vector<minimizer_t>& t_minimizers,
                       const std::unique_ptr<fastaq::FastAQ>& reference, const paired_reads_t& paired_reads,
                       mapping_params_t& parameters) {
  mapping_params_t temp = parameters;
  temp.all = false;
  temp.band = -1;

  uint32_t size = 0.05 * paired_reads.first.size();

  std::vector<uint32_t> ins;
  ins.reserve(size);

  for (uint32_t i = 0; ins.size() < size && i < paired_reads.first.size(); ++i) {
    std::vector<mapping_t> ms1 = process_single(ref_index, t_minimizers, reference, paired_reads.first[i], temp);
    std::vector<mapping_t> ms2 = process_single(ref_index, t_minimizers, reference, paired_reads.second[i], temp);
    if (ms1.size() && ms2.size() && ms1[0].mapq >= 20 && ms2[0].mapq >= 20) {
      ins.push_back(ms1[0].pos < ms2[0].pos
                    ? (ms2[0].pos + paired_reads.second[i]->sequence.size()) - ms1[0].pos
                    : (ms1[0].pos + paired_reads.first[i]->sequence.size()) - ms2[0].pos);
    }
  }

  std::sort(ins.begin(), ins.end());

  uint32_t irq, start_val, end_val, start, end, j;
  irq = 2 * (ins[ins.size() * 3 / 4] - ins[ins.size() / 4]);
  start_val = ins[ins.size() / 4] - irq;
  end_val = ins[ins.size() * 3 / 4] + irq;
  for (j = 0; j < ins.size() && ins[j] < start_val; ++j);
  start = j;
  for (j = ins.size() - 1; j > 0 && ins[j] > end_val; --j);
  end = j;
  
  uint32_t sum = 0;
  uint32_t count = 0;
  for (uint32_t i = start; i < end; ++i) {
    sum += ins[i];
    count++;
  }
  parameters.insert_size = (uint32_t)round((float)sum / count);

  sum = 0;
  for (uint32_t i = start; i < end; ++i) {
    sum += (ins[i] - parameters.insert_size) * (ins[i] - parameters.insert_size);
  }
  parameters.sd = sqrt(sum / (count - 1.0f));
}

// Map single reads
// Args: ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       reads        - list of FastAQ representations of reads
//       parameters   - mapping parameters
//       t_start      - starting index of reads to be mapped       // used for parallelization
//       t_end        - last index of reads to be mapped           // [t_start, t_end)
// Return: mapping result in SAM format
std::string map_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                       const std::vector<minimizer_t>& t_minimizers,
                       const std::unique_ptr<fastaq::FastAQ>& reference, 
                       const std::vector<std::unique_ptr<fastaq::FastAQ>>& reads, const mapping_params_t& parameters,
                       uint32_t t_start, uint32_t t_end) {
  std::string sam;

  for (uint32_t i = t_start; i < t_end; ++i) {
    std::vector<mapping_t> ms = process_single(ref_index, t_minimizers, reference, reads[i], parameters);

    if (ms.size()) {
      for (const auto& m : ms) {
        sam += sam_format(m);
      }
    } else {
      sam += unmapped_sam(reads[i]->name, reads[i]->sequence, reads[i]->quality, 0, 0, 0);
    }
  }

  return sam;
}

// Map paired-end reads as single reads
// Args: ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       paired_reads - pair of lists of FastAQ representations of paired reads
//       parameters   - mapping parameters
//       t_start      - starting index of reads to be mapped       // used for parallelization as:
//       t_end        - last index of reads to be mapped           // [t_start, t_end)
// Return: mapping result in SAM format
std::string map_as_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                          const std::vector<minimizer_t>& t_minimizers,
                          const std::unique_ptr<fastaq::FastAQ>& reference, const paired_reads_t& paired_reads,
                          const mapping_params_t& parameters, uint32_t t_start, uint32_t t_end) {
  std::string sam;

  for (uint32_t i = t_start; i < t_end; ++i) {
    std::vector<mapping_t> ms1 = process_single(ref_index, t_minimizers, reference, paired_reads.first[i], parameters);
    std::vector<mapping_t> ms2 = process_single(ref_index, t_minimizers, reference, paired_reads.second[i], parameters);

    if (ms1.size()) {
      for (const auto& m : ms1) {
        sam += sam_format(m);
      }
    } else {
      sam += unmapped_sam(paired_reads.first[i]->name, paired_reads.first[i]->sequence, paired_reads.first[i]->quality, 0, 0, 0);
    }
    if (ms2.size()) {
      for (const auto& m : ms2) {
        sam += sam_format(m);
      }
    } else {
      sam += unmapped_sam(paired_reads.second[i]->name, paired_reads.second[i]->sequence, paired_reads.second[i]->quality, 0, 0, 0);
    }
  }

  return sam;
}

// Map paired-end reads as single reads
// Args: ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       paired_reads - pair of lists of FastAQ representations of paired reads
//       parameters   - mapping parameters
//       t_start      - starting index of reads to be mapped       // used for parallelization as:
//       t_end        - last index of reads to be mapped           // [t_start, t_end)
// Return: mapping result in SAM format
std::string map_paired(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                       const std::vector<minimizer_t>& t_minimizers,
                       const std::unique_ptr<fastaq::FastAQ>& reference, const paired_reads_t& paired_reads,
                       const mapping_params_t& parameters, uint32_t t_start, uint32_t t_end) {
  std::string sam;

  for (uint32_t i = t_start; i < t_end; ++i) {
    std::pair<int32_t, int32_t> clipped1 = clip(paired_reads.first[i]->sequence, parameters.k, parameters.w);
    std::pair<int32_t, int32_t> clipped2 = clip(paired_reads.second[i]->sequence, parameters.k, parameters.w);
    if (clipped1.first < 0 || clipped2.first < 0) {
      sam += unmapped_sam(paired_reads.first[i]->name, paired_reads.first[i]->sequence,
                          paired_reads.first[i]->quality, 1, 1, 0)
             + unmapped_sam(paired_reads.second[i]->name, paired_reads.second[i]->sequence, 
                            paired_reads.second[i]->quality, 1, 0, 1);
      continue;
    }

    paired_minimizers_t p_minimizers(brown::minimizers(paired_reads.first[i]->sequence.c_str() + clipped1.first,
                                                       paired_reads.first[i]->sequence.size() - clipped1.first - clipped1.second,
                                                       parameters.k, parameters.w),
                                     brown::minimizers(paired_reads.second[i]->sequence.c_str() + clipped2.first,
                                                       paired_reads.second[i]->sequence.size() - clipped2.first - clipped2.second,
                                                       parameters.k, parameters.w));
    if (clipped1.first) {
      for (auto& min : p_minimizers.first) {
        std::get<1>(min) += clipped1.first;
      }
    }
    if (clipped2.first) {
      for (auto& min : p_minimizers.second) {
        std::get<1>(min) += clipped2.first;
      }
    }

    split_hits_t hits1;
    split_hits_t hits2;
    find_minimizer_hits(hits1.first, hits1.second, ref_index, t_minimizers, p_minimizers.first);
    find_minimizer_hits(hits2.first, hits2.second, ref_index, t_minimizers, p_minimizers.second);

    radixsort(hits1.first);
    radixsort(hits1.second);
    radixsort(hits2.first);
    radixsort(hits2.second);

    std::pair<bin_t, bin_t> candidates1(
        extract_candidates(hits1.first, parameters.threshold, paired_reads.first[i]->sequence.size()),
        extract_candidates(hits2.second, parameters.threshold, paired_reads.second[i]->sequence.size()));
    std::pair<bin_t, bin_t> candidates2(
        extract_candidates(hits1.second, parameters.threshold, paired_reads.first[i]->sequence.size()),
        extract_candidates(hits2.first, parameters.threshold, paired_reads.second[i]->sequence.size()));

    paired_checked_t checked1 = check_pairing(candidates1, parameters.insert_size,
                                              paired_reads.second[i]->sequence.size());
    paired_checked_t checked2 = check_pairing(candidates2, parameters.insert_size,
                                              paired_reads.first[i]->sequence.size());

    std::vector<std::pair<mapping_t, mapping_t>> mappings;
    process_pairs(mappings, checked1, reference, paired_reads.first[i], paired_reads.second[i], parameters, clipped1, clipped2);
    process_pairs(mappings, checked2, reference, paired_reads.first[i], paired_reads.second[i], parameters, clipped1, clipped2);

    if (mappings.size() == 0) {
      sam += unmapped_sam(paired_reads.first[i]->name, paired_reads.first[i]->sequence,
                          paired_reads.first[i]->quality, 1, 1, 0)
             + unmapped_sam(paired_reads.second[i]->name, paired_reads.second[i]->sequence, 
                            paired_reads.second[i]->quality, 1, 0, 1);
      continue;
    }

    std::sort(mappings.begin(), mappings.end(), 
              [] (const std::pair<mapping_t, mapping_t>& a, 
                  const std::pair<mapping_t, mapping_t>& b) {
                return a.first.mapq + a.second.mapq > b.first.mapq + b.second.mapq;
              }
    );

    std::string sam1;
    std::string sam2;

    if (parameters.all) {
      for (uint32_t j = 0; j < mappings.size(); ++j) {
        mappings[j].first.mapq /= mappings.size();
        mappings[j].second.mapq /= mappings.size();

        if (j) {
          mappings[j].first.flag |= 0x100;
          mappings[j].second.flag |= 0x100;
          mappings[j].first.mapq = 0;
          mappings[j].second.mapq = 0;
        }
        
        sam1 += sam_format(mappings[j].first);
        sam2 += sam_format(mappings[j].second);
      }
    } else {
      uint32_t j = 1;
      int32_t as = mappings[0].first.as + mappings[0].second.as;

      while (j < mappings.size() && mappings[j].first.as + mappings[j].second.as == as) ++j;
      uint32_t pick = std::rand() % j;
      
      mappings[pick].first.mapq /= mappings.size();
      mappings[pick].second.mapq /= mappings.size();

      sam1 = sam_format(mappings[pick].first);
      sam2 = sam_format(mappings[pick].second);
    }

    sam += sam1 + sam2;
  }

  return sam;
}