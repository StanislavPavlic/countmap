#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <algorithm>
#include <iostream>

#include "sam.hpp"
#include "mapping_params.hpp"
#include "mapping.hpp"
#include "ksw2.h"

std::unordered_map<uint8_t, uint8_t> c = {{'C', 0}, {'A', 1}, {'T', 2}, {'U', 2}, {'G', 3}};

// Update offset values according to clipping and strand
// Args: clip   - clipping size
//       off1   - first offset value
//       off2   - second offset value
//       strand - mapping strand
// Return: none
void update_offset(uint32_t clip, uint32_t& off1, uint32_t& off2, bool strand) {
  if (clip) {
    if (strand) {
      if (clip > off1) {
        off1 = clip;
      }
    } else {
      if (clip > off2) {
        off2 = clip;
      }
    }
  }
}

// Perform KSW2 algorithm on found region, return cigar
// Args: target     - target sequence
//       query      - query sequence
//       region     - region on query and target
//       parameters - wrapper of significant mapping parameters
// Return: number of exact matches, alignment score, cigar string
std::tuple<uint32_t, int32_t, std::string> ksw2(const char* target, const uint32_t t_len, 
                                                const char* query, const uint32_t q_len, 
                                                const mapping_params_t& parameters) {
  int a = parameters.mch > 0 ? parameters.mch : -parameters.mch;
  int b = parameters.mis < 0 ? parameters.mis : -parameters.mis;
  int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
  uint8_t *ts, *qs;
  ksw_extz_t ez;

  memset(&ez, 0, sizeof(ksw_extz_t));
  ts = (uint8_t*)malloc(t_len);
  qs = (uint8_t*)malloc(q_len);

  c['N'] = std::rand() % 4;
  for (uint32_t i = 0; i < t_len; ++i) {
    ts[i] = c[(uint8_t)target[i]];
  }
  for (uint32_t i = 0; i < q_len; ++i) {
    qs[i] = c[(uint8_t)query[i]];
  }

  ksw_extz2_sse(0, q_len, qs, t_len, ts, 5, mat, 
                parameters.gapo, parameters.gape, 
                (parameters.band == -1 ? 0.1 * q_len : parameters.band), -1, 0, 0, &ez); // Set band width to 10% of query length

  std::string cigar;
  uint32_t matches = 0;
  uint32_t t_pos = 0;
  uint32_t q_pos = 0;
  for (int i = 0; i < ez.n_cigar; ++i) {
    if ("MIDNSHP=X"[ez.cigar[i]&0xf] == 'M') {
      for (uint32_t j = 0; j < (ez.cigar[i]>>4); ++j) {
        if (target[t_pos + j] == query[q_pos + j]) {
          matches++;
        }
      }
      t_pos += ez.cigar[i]>>4;
      q_pos += ez.cigar[i]>>4;
    } else if ("MIDNSHP=X"[ez.cigar[i]&0xf] == 'I') {
      q_pos += ez.cigar[i]>>4;
    } else if ("MIDNSHP=X"[ez.cigar[i]&0xf] == 'D') {
      t_pos += ez.cigar[i]>>4;
    }

    cigar += std::to_string(ez.cigar[i]>>4) + "MIDNSHP=X"[ez.cigar[i]&0xf];
  }

  free(ez.cigar); free(ts); free(qs);

  return std::make_tuple(matches, ez.score, cigar);
}

// Convert mapping to SAM formatted string
// Args: m - mapping information
// Return: SAM format string
std::string sam_format(const mapping_t& m) {
  std::string sam = m.qname + "\t" +
                    std::to_string(m.flag) + "\t" +
                    m.rname + "\t" +
                    std::to_string(m.pos) + "\t" +
                    std::to_string(m.mapq) + "\t" +
                    m.cigar + "\t" +
                    m.rnext + "\t" +
                    std::to_string(m.pnext) + "\t" +
                    std::to_string(m.tlen) + "\t" +
                    (m.flag & 0x100 ? "*" : m.seq) + "\t" +
                    (m.flag & 0x100 ? "*" : m.qual) + "\t" +
                    "NM:i:" + std::to_string(m.nm) + "\t" +
                    "AS:i:" + std::to_string(m.as) + "\n";
  return sam;
}

// Generates mapping in SAM format
// Args: qname      - query name from read file
//       query      - query sequence
//       qual       - query sequence per base qualities
//       rname      - reference name from reference file
//       ref        - reference sequence
//       region     - mapping region
//       parameters - mapping parameters
//       clipped    - query sequence clipping
// Return: mapping
mapping_t single_mapping(const std::string& qname, const std::string& query, 
                         const std::string& qual, const std::string& rname, 
                         const std::string& ref, const region_t& region, 
                         const mapping_params_t& parameters, 
                         const std::pair<int32_t, int32_t>& clipped) {
  uint32_t start_off = std::get<2>(region.first)
                        ? query.size() - std::get<0>(region.second)
                        : std::get<0>(region.first);
  uint32_t end_off = std::get<2>(region.first)
                      ? std::get<0>(region.first)
                      : query.size() - std::get<0>(region.second);
  
  std::string preclip, postclip;
  update_offset(clipped.first, end_off, start_off, std::get<2>(region.first));
  update_offset(clipped.second, start_off, end_off, std::get<2>(region.first));
  
  uint32_t len = (query.size() - end_off) - start_off;

  std::tuple<uint32_t, int32_t, std::string> cigar = ksw2(ref.c_str() + std::get<1>(region.first) + start_off, len, 
                                                          query.c_str() + start_off, len, 
                                                          parameters);

  mapping_t m;

  m.qname = qname.substr(0, qname.find('/', 0));
  m.flag = std::get<2>(region.first) ? 0x10 : 0x0;
  m.rname = rname;
  m.pos = std::get<1>(region.first) + 1;
  m.mapq = (uint32_t)round((double)std::get<0>(cigar) / (query.size() - start_off - end_off) * 60);;
  m.cigar = preclip + std::get<2>(cigar) + postclip;
  m.rnext = "*";
  m.pnext = 0;
  m.tlen = 0;
  m.seq = query;
  m.qual = qual;
  m.nm = query.size() - std::get<0>(cigar) - start_off - end_off;
  m.as = std::get<1>(cigar);

  return m;
}

// Generates paired mapping in SAM format
// Args: qname       - query name from read file
//       query1      - query sequence from first read
//       qual1       - query sequence per base qualities from first read
//       query2      - query sequence from second read
//       qual2       - query sequence per base qualities from second read
//       rname       - reference name from reference file
//       ref         - reference sequence
//       region_pair - paired mapping region
//       parameters  - mapping parameters
//       clipped1    - query1 sequence clipping
//       clipped2    - query2 sequence clipping
// Return: paired mapping
std::pair<mapping_t, mapping_t> pair_mapping(const std::string& qname,
    const std::string& query1, const std::string& qual1, const std::string& query2, const std::string& qual2,
    const std::string& rname, const std::string& ref,
    const std::pair<region_t, region_t>& region_pair, const mapping_params_t& parameters,
    const std::pair<int32_t, int32_t>& clipped1, 
    const std::pair<int32_t, int32_t>& clipped2) {
  int32_t insert_size = std::get<1>(region_pair.first.first) < std::get<1>(region_pair.second.first)
                        ? std::get<1>(region_pair.second.second) - std::get<1>(region_pair.first.first)
                        : std::get<1>(region_pair.second.first) - std::get<1>(region_pair.first.second);

  uint32_t start_off1 = std::get<2>(region_pair.first.first)
                        ? query1.size() - std::get<0>(region_pair.first.second)
                        : std::get<0>(region_pair.first.first);
  uint32_t start_off2 = std::get<2>(region_pair.second.first)
                        ? query2.size() - std::get<0>(region_pair.second.second)
                        : std::get<0>(region_pair.second.first);
  uint32_t end_off1 = std::get<2>(region_pair.first.first)
                      ? std::get<0>(region_pair.first.first)
                      : query1.size() - std::get<0>(region_pair.first.second);
  uint32_t end_off2 = std::get<2>(region_pair.second.first)
                      ? std::get<0>(region_pair.second.first)
                      : query2.size() - std::get<0>(region_pair.second.second);
  update_offset(clipped1.first, end_off1, start_off1, std::get<2>(region_pair.first.first));
  update_offset(clipped1.second, start_off1, end_off1, std::get<2>(region_pair.first.first));
  update_offset(clipped2.first, end_off2, start_off2, std::get<2>(region_pair.second.first));
  update_offset(clipped2.second, start_off2, end_off2, std::get<2>(region_pair.second.first));
  
  std::string preclip1, postclip1, preclip2, postclip2;
  if (start_off1) {
    preclip1 = std::to_string(start_off1) + "S";
  }
  if (end_off1) {
    postclip1 = std::to_string(end_off1) + "S";
  }
  if (start_off2) {
    preclip2 = std::to_string(start_off2) + "S";
  }
  if (end_off2) {
    postclip2 = std::to_string(end_off2) + "S";
  }
  
  uint32_t len1 = (query1.size() - end_off1) - start_off1;
  uint32_t len2 = (query2.size() - end_off2) - start_off2;

  std::tuple<uint32_t, int32_t, std::string> cigar1 = ksw2(ref.c_str() + std::get<1>(region_pair.first.first) + start_off1, len1, 
                                                           query1.c_str() + start_off1, len1, 
                                                           parameters);
  std::tuple<uint32_t, int32_t, std::string> cigar2 = ksw2(ref.c_str() + std::get<1>(region_pair.second.first) + start_off2, len2, 
                                                           query2.c_str() + start_off2, len2, 
                                                           parameters);

  int prop_aligned = (float)abs(abs(insert_size) - (int32_t)parameters.insert_size) < 5 * parameters.sd ? 0x2 : 0x0;

  double z = ((double)abs(insert_size) - parameters.insert_size) / (5 * parameters.sd);              

  mapping_t m1, m2;

  m1.qname = m2.qname = qname.substr(0, qname.find('/', 0));
  m1.rname = m2.rname = rname;
  m1.rnext = m2.rnext = "=";

  m1.flag = 0x1 | prop_aligned | (std::get<2>(region_pair.first.first) ? 0x10 : 0x0) 
                | (std::get<2>(region_pair.second.first) ? 0x20 : 0x0) | 0x40;
  m1.pos = std::get<1>(region_pair.first.first) + 1;
  m1.mapq = (uint32_t)round(std::max((double)std::get<0>(cigar1) / (query1.size() - start_off1 - end_off1) * 60 - z*z, (double)0));
  m1.cigar = preclip1 + std::get<2>(cigar1) + postclip1;
  m1.pnext = std::get<1>(region_pair.second.first) + 1;
  m1.tlen = insert_size;
  m1.seq = query1;
  m1.qual = qual1;
  m1.nm = query1.size() - std::get<0>(cigar1) - start_off1 - end_off1;
  m1.as = std::get<1>(cigar1);
  
  m2.flag = 0x1 | prop_aligned | (std::get<2>(region_pair.second.first) ? 0x10 : 0x0) 
                | (std::get<2>(region_pair.first.first) ? 0x20 : 0x0) | 0x80;
  m2.pos = std::get<1>(region_pair.second.first) + 1;
  m2.mapq = (uint32_t)round(std::max((double)std::get<0>(cigar2) / (query2.size() - (start_off2 + end_off2)) * 60 - z*z, (double)0));
  m2.cigar = preclip2 + std::get<2>(cigar2) + postclip2;
  m2.pnext = std::get<1>(region_pair.first.first) + 1;
  m2.tlen = -insert_size;
  // if (std::get<0>(region_pair.second.first) > 0 || std::get<0>(region_pair.second.second) < query2.size()) {
  //   m2.seq = query2.substr(std::get<0>(region_pair.second.first), std::get<0>(region_pair.second.second) - std::get<0>(region_pair.second.first));
  //   m2.qual = qual2.substr(std::get<0>(region_pair.second.first), std::get<0>(region_pair.second.second) - std::get<0>(region_pair.second.first));
  // } else {
    m2.seq = query2;
    m2.qual = qual2;
  // }
  m2.nm = query2.size() - std::get<0>(cigar2) - start_off2 - end_off2;
  m2.as = std::get<1>(cigar2);

  return std::make_pair(m1, m2);
}

// Generates SAM format output for unmapped read
// Args: qname - query name from read file
//       query - query sequence
//       qual  - query sequence per base qualities
//       pair  - paired reads flag
//       first - first in segment flag
//       last  - last in segment flag
// Return: unmapped read SAM format output
std::string unmapped_sam(const std::string& qname, const std::string& query, const std::string& qual, 
                         bool pair, bool first, bool last) {
  std::string sam_name = qname.substr(0, qname.find('/', 0));

  int flag = 0x4;
  if (pair) {
    flag |= 0x1 | 0x8;
    if (first) flag |= 0x40;
    if (last) flag |= 0x80;
  }

  std::string sam = sam_name + "\t" +
                    std::to_string(flag) + "\t" +
                    "*" + "\t" +
                    "0" + "\t" +
                    "255" + "\t" +
                    "*" + "\t" +
                    "*" + "\t" +
                    "0" + "\t" +
                    "0" + "\t" +
                    query + "\t" +
                    qual + "\n";

  return sam;
}