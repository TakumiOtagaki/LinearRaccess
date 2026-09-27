#include "beam_inside_outside.hpp"
#include "linearraccess/dp_table_api.hpp"

#include <algorithm>
#include <cmath>
vector<double> LinCapR::calc_accessibility(int window) const {
	return calc_accessibility_by_loop(window).total;
}

vector<vector<double>> LinCapR::calc_accessibility(const vector<int>& windows) const {
	auto by = calc_accessibility_by_loop(windows);
	vector<vector<double>> totals;
	totals.reserve(by.size());
	for(auto& b : by) totals.push_back(std::move(b.total));
	return totals;
}

AccessibilityByLoop LinCapR::calc_accessibility_by_loop(int window) const {
	if(window <= 0 || window > seq_n) return {};
	return calc_accessibility_by_loop(std::vector<int>{window}).front();
}

// All windows in one pass: every loop is enumerated (and its energy and
// Boltzmann weight evaluated) once, and its mass is added to each window's
// difference array.  Per window, the additions happen in the same order as a
// single-window pass, so the result is bit-identical to calling the
// single-window version once per window.
std::vector<AccessibilityByLoop>
LinCapR::calc_accessibility_by_loop(const std::vector<int>& requested) const {
	// Windows outside [1, seq_n] get an empty result, as in the
	// single-window version.
	std::vector<AccessibilityByLoop> result(requested.size());
	std::vector<int> windows, slot;
	for(size_t r = 0; r < requested.size(); r++){
		if(requested[r] > 0 && requested[r] <= seq_n){
			windows.push_back(requested[r]);
			slot.push_back(static_cast<int>(r));
		}
	}
	const int nw = static_cast<int>(windows.size());
	if(nw == 0) return result;
	std::vector<AccessibilityByLoop> outs(nw);

	const double logZ = alpha_O[seq_n - 1];
	int min_window = windows[0];
	for(const int window : windows) min_window = std::min(min_window, window);

	std::vector<std::vector<double>> diff_hairpin(nw), diff_bulge(nw), diff_internal(nw), diff_multi(nw);
	for(int w = 0; w < nw; w++){
		const int starts = seq_n - windows[w] + 1;
		AccessibilityByLoop& out = outs[w];
		out.total.assign(starts, 0.0);
		out.exterior.assign(starts, 0.0);
		out.hairpin.assign(starts, 0.0);
		out.bulge.assign(starts, 0.0);
		out.internal.assign(starts, 0.0);
		out.multiloop.assign(starts, 0.0);
		diff_hairpin[w].assign(starts + 1, 0.0);
		diff_bulge[w].assign(starts + 1, 0.0);
		diff_internal[w].assign(starts + 1, 0.0);
		diff_multi[w].assign(starts + 1, 0.0);
	}

	// Adds mass to every window of length window[w] inside [s, e], for each w
	// with window[w] <= e - s + 1.
	auto range_add = [&](std::vector<std::vector<double>>& diffs, int s, int e, double mass) {
		const int seg_len = e - s + 1;
		for(int w = 0; w < nw; w++){
			const int window = windows[w];
			if(seg_len < window) continue;
			const int max_start = seq_n - window;
			int l = s;
			int r = e - window + 1;
			if(l < 0) l = 0;
			if(r > max_start) r = max_start;
			if(l <= r){
				diffs[w][l] += mass;
				diffs[w][r + 1] -= mass;
			}
		}
	};

	// (A) Exterior contribution using alpha_O / beta_O and external unpaired energies.
	std::vector<double> ext_prefix(seq_n + 1, 0.0);
	for(int i = 0; i < seq_n; i++){
		ext_prefix[i + 1] = ext_prefix[i] + _energy->energy_external_unpaired(i, i) / _energy->kT();
	}
	for(int w = 0; w < nw; w++){
		const int window = windows[w];
		for(int i = 0; i <= seq_n - window; i++){
			const double logW = (i > 0 ? alpha_O[i - 1] : 0.0)
				+ (i + window < seq_n ? beta_O[i + window] : 0.0)
				- (ext_prefix[i + window] - ext_prefix[i]);
			outs[w].exterior[i] = std::exp(logW - logZ);
		}
	}

	// (B) Hairpin segments from beta_SE.
	for(int e = 0; e < seq_n; e++){
		for(const auto& kv : beta_SE[e]){
			const int s = kv.first;
			const double score = kv.second;
			if(e - s + 1 > c_hairpin) continue;
			if(e - s + 1 < min_window) continue;
			const int outer_i = s - 1;
			const int outer_j = e + 1;
			const double logW = score - _energy->energy_hairpin(outer_i, outer_j) / _energy->kT();
			const double mass = std::exp(logW - logZ);
			range_add(diff_hairpin, s, e, mass);
		}
	}

	// (C) Internal/Bulge segments via beta_SE + alpha_S.  A loop contributes
	// only through an unpaired side of length >= some window, so loops whose
	// both sides are shorter than the smallest window are skipped before the
	// energy is evaluated.
	for(int e = 0; e < seq_n; e++){
		for(const auto& kv : beta_SE[e]){
			const int s = kv.first;
			const double score = kv.second;
			for(int p = s; p <= std::min(s + MAXLOOP, e - 1); p++){
				for(int q = e; q >= p + TURN + 1 && (p - s) + (e - q) <= MAXLOOP; q--){
					if((p == s && q == e) || !lcr::dp::contains(alpha_S, p, q)) continue;
					const bool left = p - 1 >= s && (p - s) >= min_window;
					const bool right = q + 1 <= e && (e - q) >= min_window;
					if(!left && !right) continue;
					auto it_a = alpha_S[q].find(p);
					if(it_a == alpha_S[q].end()) continue;
					const double logW = score + it_a->second
						- _energy->energy_loop(s - 1, e + 1, p, q) / _energy->kT();
					const double mass = std::exp(logW - logZ);

					// Left unpaired segment [s, p-1].
					if(left){
						range_add(q == e ? diff_bulge : diff_internal, s, p - 1, mass);
					}

					// Right unpaired segment [q+1, e].
					if(right){
						range_add(p == s ? diff_bulge : diff_internal, q + 1, e, mass);
					}
				}
			}
		}
	}

	// (D1) Multiloop left-unpaired via beta_M + alpha_MB.
	for(int k = 0; k < seq_n; k++){
		for(const auto& kv : alpha_MB[k]){
			const int p = kv.first;
			const double score = kv.second;
			for(int j = p - min_window; j >= std::max(0, p - c_multi); j--){
				if(!lcr::dp::contains(beta_M, j, k)) continue;
				auto it_b = beta_M[k].find(j);
				if(it_b == beta_M[k].end()) continue;
				const double logW = score + it_b->second
					- _energy->energy_multi_unpaired(j, p - 1) / _energy->kT();
				const double mass = std::exp(logW - logZ);
				range_add(diff_multi, j, p - 1, mass);
			}
		}
	}

	// (D2) Multiloop right-unpaired via beta_M2 + alpha_S.
	for(int q = 0; q < seq_n; q++){
		for(const auto& kv : alpha_S[q]){
			const int j = kv.first;
			const double score = kv.second;
			for(int k = q + min_window; k <= std::min(seq_n - 1, q + c_multi); k++){
				if(!lcr::dp::contains(beta_M2, j, k)) continue;
				auto it_b = beta_M2[k].find(j);
				if(it_b == beta_M2[k].end()) continue;
				const double logW = score + it_b->second
					- (_energy->energy_multi_bif(j, q) + _energy->energy_multi_unpaired(q + 1, k)) / _energy->kT();
				const double mass = std::exp(logW - logZ);
				range_add(diff_multi, q + 1, k, mass);
			}
		}
	}

	for(int w = 0; w < nw; w++){
		AccessibilityByLoop& out = outs[w];
		double running_hairpin = 0.0;
		double running_bulge = 0.0;
		double running_internal = 0.0;
		double running_multi = 0.0;
		for(int i = 0; i <= seq_n - windows[w]; i++){
			running_hairpin += diff_hairpin[w][i];
			running_bulge += diff_bulge[w][i];
			running_internal += diff_internal[w][i];
			running_multi += diff_multi[w][i];

			// Preserve raw posterior sums.  Silently clamping numerical violations
			// would hide approximation error from validation and downstream users.
			out.hairpin[i] = running_hairpin;
			out.bulge[i] = running_bulge;
			out.internal[i] = running_internal;
			out.multiloop[i] = running_multi;
			out.total[i] = out.exterior[i] + out.hairpin[i] + out.bulge[i] + out.internal[i] + out.multiloop[i];
		}
	}

	for(int w = 0; w < nw; w++) result[slot[w]] = std::move(outs[w]);
	return result;
}
