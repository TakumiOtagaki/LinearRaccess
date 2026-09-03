#pragma once

#include "linearraccess/miscs.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace lcr {
namespace io {

void loop_profile_io(std::ostream& os,
                     const std::string& seq_name,
                     const std::vector<Float>& prob_B,
                     const std::vector<Float>& prob_E,
                     const std::vector<Float>& prob_H,
                     const std::vector<Float>& prob_I,
                     const std::vector<Float>& prob_M,
                     const std::vector<Float>& prob_S);

void accessibility_raccess_io(std::ostream& os,
                              const std::string& seq_name,
                              int seq_len,
                              const std::vector<int>& access_lens,
                              const std::vector<std::vector<double>>& probs,
                              double rt);

void accessibility_probability_io(std::ostream& os,
                                   const std::string& seq_name,
                                   int seq_len,
                                   const std::vector<int>& access_lens,
                                   const std::vector<std::vector<double>>& probs);

void accessibility_raccess_byloop_io(std::ostream& os,
                                     const std::string& seq_name,
                                     int seq_len,
                                     const std::vector<int>& access_lens,
                                     const std::vector<std::vector<double>>& total_probs,
                                     const std::vector<std::vector<double>>& exterior_probs,
                                     const std::vector<std::vector<double>>& hairpin_probs,
                                     const std::vector<std::vector<double>>& bulge_probs,
                                     const std::vector<std::vector<double>>& internal_probs,
                                     const std::vector<std::vector<double>>& multiloop_probs,
                                     double rt);

} // namespace io
} // namespace lcr
