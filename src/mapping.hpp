#pragma once

#include <cstdint>
#include <string>

struct mapping_t {
  std::string qname;
  int flag;
  std::string rname;
  uint32_t pos;
  uint32_t mapq;
  std::string cigar;
  std::string rnext;
  uint32_t pnext;
  int32_t tlen;
  std::string seq;
  std::string qual;
  uint32_t nm;
  int32_t as;
};