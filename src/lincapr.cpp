#include "beam_inside_outside.hpp"
#include "linearraccess/dp_table_api.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

// calc structural profile
void LinCapR::calc_profile(){
	if(debug_outer_alpha_log){
		int finite = 0;
		for(int i = 0; i < seq_n; i++){
			if(alpha_O[i] > -INF / 2) finite++;
		}
		cerr << "debug_outer_alpha_state: phase=before_profile"
		     << " alpha_O0=" << (seq_n > 0 ? alpha_O[0] : 0.0)
		     << " alpha_O_last=" << (seq_n > 0 ? alpha_O[seq_n - 1] : 0.0)
		     << " finite=" << finite
		     << endl;
	}
	const Float logZ = alpha_O[seq_n - 1];

	for(int k = 0; k < seq_n; k++){
		for(const auto [j, score] : beta_SE[k]){
			// H: unpaired range is derived from raw closing pair (outer_i, outer_j).
			const int outer_i = j - 1;
			const int outer_j = k + 1;
			const int unp_l = outer_i + 1;
			const int unp_r = outer_j - 1;
			lcr::dp::add_range(prob_H, unp_l, unp_r,
			                   exp(score - _energy->energy_hairpin(outer_i, outer_j) / _energy->kT() - logZ));

			// B, I
			for(int p = j; p <= min(j + MAXLOOP, k - 1); p++){
				for(int q = k; q >= p + TURN + 1 && (p - j) + (k - q) <= MAXLOOP; q--){
					if((p == j && q == k) || !lcr::dp::contains(alpha_S, p, q)) continue;
					auto it_a = alpha_S[q].find(p);
					if(it_a == alpha_S[q].end()) continue;
					const Float new_score = exp(score + it_a->second - _energy->energy_loop(j - 1, k + 1, p, q) / _energy->kT() - logZ);
					lcr::dp::add_range((q == k ? prob_B : prob_I), j, p - 1, new_score);
					lcr::dp::add_range((p == j ? prob_B : prob_I), q + 1, k, new_score);
				}
			}
		}
	}
	lcr::dp::prefix_sum(prob_B);
	lcr::dp::prefix_sum(prob_H);
	lcr::dp::prefix_sum(prob_I);

	// M
	for(int k = 0; k < seq_n; k++){
		for(const auto [p, score] : alpha_MB[k]){
			for(int j = p - 1; j >= max(0, p - c_multi); j--){
				if(!lcr::dp::contains(beta_M, j, k)) continue;
				const Float new_score = exp(score + beta_M[k][j] - _energy->energy_multi_unpaired(j, p - 1) / _energy->kT() - logZ);
				lcr::dp::add_range(prob_M, j, p - 1, new_score);
			}
		}
	}
	for(int q = 0; q < seq_n; q++){
		for(const auto [j, score] : alpha_S[q]){
			for(int k = q + 1; k <= min(seq_n - 1, q + c_multi); k++){
				if(!lcr::dp::contains(beta_M2, j, k)) continue;
				const Float new_score = exp(score + beta_M2[k][j] - (_energy->energy_multi_bif(j, q) + _energy->energy_multi_unpaired(q + 1, k)) / _energy->kT() - logZ);
				lcr::dp::add_range(prob_M, q + 1, k, new_score);
			}
		}
	}
	lcr::dp::prefix_sum(prob_M);

	// S
	for(int j = 0; j < seq_n; j++){
		for(const auto [i, score] : alpha_S[j]){
			auto it_b = beta_S[j].find(i);
			if(it_b == beta_S[j].end()) continue;
			const Float new_score = exp(score + it_b->second - logZ);
			prob_S[i] += new_score;
			prob_S[j] += new_score;
		}
	}

	// E
	for(int i = 0; i < seq_n; i++){
		const double prefix = (i > 0 ? alpha_O[i - 1] : 0.0);
		const double suffix = (i + 1 < seq_n ? beta_O[i + 1] : 0.0);
		prob_E[i] = exp(prefix + suffix
		                - _energy->energy_external_unpaired(i, i) / _energy->kT()
		                - logZ);
	}

	// regularize
	for(int i = 0; i < seq_n; i++){
		Float sum_prob_i = 0;
		for(int j = 0; j < NPROBS; j++){
			// negative probabilities to 0
			if(probs[j]->at(i) < 0) probs[j]->at(i) = 0;
			sum_prob_i += probs[j]->at(i);
		}
		if(fabs(sum_prob_i - 1.0) > normalize_warn_eps){
			cerr << "warn: prob_sum[" << i << "]=" << sum_prob_i << endl;
		}
		if(!normalize_profiles) continue;
		// sum of probabilities to 1
		for(int j = 0; j < NPROBS; j++) probs[j]->at(i) /= sum_prob_i;
	}
	if(debug_outer_alpha_log){
		int finite = 0;
		for(int i = 0; i < seq_n; i++){
			if(alpha_O[i] > -INF / 2) finite++;
		}
		cerr << "debug_outer_alpha_state: phase=after_profile"
		     << " alpha_O0=" << (seq_n > 0 ? alpha_O[0] : 0.0)
		     << " alpha_O_last=" << (seq_n > 0 ? alpha_O[seq_n - 1] : 0.0)
		     << " finite=" << finite
		     << endl;
	}
}
