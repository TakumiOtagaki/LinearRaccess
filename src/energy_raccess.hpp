/*
 * Raccess-backed energy model for LinearCapR-style DP.
 */
#pragma once

#include "linearraccess/energy_api.hpp"
#include "energy_model.hpp"
#include <memory>

namespace lcr {

std::unique_ptr<EnergyApi> make_raccess_energy(const energy::Params* params = nullptr,
                                               bool use_linearcapr_external = false);
double raccess_rt_kcal_mol();
double raccess_length_factor();
std::vector<std::vector<double>> compute_raccess_accessibility(
    const std::string& seq, const std::vector<int>& access_lens, int max_span,
    double length_factor = 0.0);
std::vector<double> compute_raccess_unpaired_1(const std::string& seq, int max_span);
void debug_raccess_multi_unpaired(const std::string& seq, int i, int j);
void debug_raccess_local(const std::string& seq, int i, int j, bool has_loop, int p, int q);
void debug_raccess_outer_branch_terms(const std::string& seq, int i, int j);
double compute_raccess_logz(const std::string& seq, int max_span);

} // namespace lcr
