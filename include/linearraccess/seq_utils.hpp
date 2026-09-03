/*
 * Sequence encoding and pairing utilities (LinearCapR style).
 */
#pragma once

#include "linearraccess/miscs.hpp"

#include <stdexcept>
#include <string>

namespace lcr {
namespace seq {

inline int base_to_num(const char base) { return ::base_to_num(base); }

inline char canonical_base(const char base) {
  switch(base) {
  case 'A': case 'a': return 'A';
  case 'C': case 'c': return 'C';
  case 'G': case 'g': return 'G';
  case 'T': case 't': case 'U': case 'u': return 'U';
  default: return '\0';
  }
}

inline std::string normalize_sequence(const std::string& sequence) {
  if(sequence.empty()) {
    throw std::invalid_argument("sequence must not be empty");
  }
  std::string normalized;
  normalized.reserve(sequence.size());
  for(std::size_t i = 0; i < sequence.size(); ++i) {
    const char base = canonical_base(sequence[i]);
    if(base == '\0') {
      throw std::invalid_argument(
          "unsupported nucleotide at 1-based position " + std::to_string(i + 1)
          + "; accepted alphabet is A,C,G,U,T (case-insensitive)");
    }
    normalized.push_back(base);
  }
  return normalized;
}

inline bool can_pair(const int a, const int b) { return (BP_pair[a][b] > 0); }

inline void build_next_pair(const vector<int>& seq_int, vector<int> next_pair[NBASE]) {
  const int seq_n = (int)seq_int.size();
  for (int i = 0; i < NBASE; ++i) next_pair[i].assign(seq_n + 1, seq_n);
  for (int i = seq_n - 1; i >= 0; --i) {
    for (int j = 0; j < NBASE; ++j) {
      next_pair[j][i] = next_pair[j][i + 1];
      if (BP_pair[seq_int[i]][j] > 0) next_pair[j][i] = i;
    }
  }
}

} // namespace seq
} // namespace lcr
