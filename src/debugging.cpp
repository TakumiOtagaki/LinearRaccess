#include "beam_inside_outside.hpp"
#ifdef LINEARRACCESS_WITH_RACCESS
#include "energy_raccess.hpp"
#endif
#include "linearraccess/dp_table_api.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <iostream>

void LinCapR::debug_stem_pairs(int idx, int topn) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_stem_pairs: idx out of range: " << idx << endl;
		return;
	}
	if(topn <= 0) topn = 1;

	struct Item {
		int i;
		int j;
		double prob;
	};
	vector<Item> items;
	items.reserve(seq_n);

	const double logZ = alpha_O[seq_n - 1];
	double total = 0.0;
	int missing_beta = 0;

	// idx as left endpoint (i = idx, j varies)
	for(int j = 0; j < seq_n; j++){
		auto it_a = alpha_S[j].find(idx);
		if(it_a == alpha_S[j].end()) continue;
		auto it_b = beta_S[j].find(idx);
		if(it_b == beta_S[j].end()){
			missing_beta++;
			continue;
		}
		const double prob = exp(it_a->second + it_b->second - logZ);
		if(prob <= 0.0) continue;
		items.push_back({idx, j, prob});
		total += prob;
	}

	// idx as right endpoint (j = idx, i varies)
	for(const auto& kv : alpha_S[idx]){
		const int i = kv.first;
		auto it_b = beta_S[idx].find(i);
		if(it_b == beta_S[idx].end()){
			missing_beta++;
			continue;
		}
		const double prob = exp(kv.second + it_b->second - logZ);
		if(prob <= 0.0) continue;
		items.push_back({i, idx, prob});
		total += prob;
	}

	sort(items.begin(), items.end(), [](const Item& a, const Item& b){
		return a.prob > b.prob;
	});

	cerr << "debug_stem_pairs i=" << idx
	     << " total=" << total
	     << " prob_S=" << prob_S[idx]
	     << " missing_beta=" << missing_beta << endl;
	const int limit = min(topn, static_cast<int>(items.size()));
	for(int k = 0; k < limit; k++){
		cerr << "  pair (" << items[k].i << "," << items[k].j << ") prob=" << items[k].prob << endl;
	}
}

void LinCapR::debug_pair(int i, int j) const {
	if(i < 0 || j < 0 || i >= seq_n || j >= seq_n || i >= j){
		cerr << "debug_pair: invalid i,j: " << i << "," << j << endl;
		return;
	}
	const double logZ = alpha_O[seq_n - 1];
	const bool pairable = can_pair(i, j);
	cerr << "debug_pair (" << i << "," << j << ")"
	     << " bases=" << seq[i] << "," << seq[j]
	     << " can_pair=" << pairable << endl;

	auto it_a = alpha_S[j].find(i);
	auto it_b = beta_S[j].find(i);
	if(it_a == alpha_S[j].end()){
		cerr << "  alpha_S missing" << endl;
	} else {
		cerr << "  alpha_S=" << it_a->second << endl;
	}
	if(it_b == beta_S[j].end()){
		cerr << "  beta_S missing" << endl;
	} else {
		cerr << "  beta_S=" << it_b->second << endl;
	}

	auto it_se_a = alpha_SE[j].find(i);
	auto it_se_b = beta_SE[j].find(i);
	if(it_se_a == alpha_SE[j].end()){
		cerr << "  alpha_SE missing" << endl;
	} else {
		cerr << "  alpha_SE=" << it_se_a->second << endl;
	}
	if(it_se_b == beta_SE[j].end()){
		cerr << "  beta_SE missing" << endl;
	} else {
		cerr << "  beta_SE=" << it_se_b->second << endl;
	}
	if(it_b != beta_S[j].end()){
		cerr << "  beta_S_norm=" << (it_b->second - logZ) << endl;
	}
	if(it_se_b != beta_SE[j].end()){
		cerr << "  beta_SE_norm=" << (it_se_b->second - logZ) << endl;
	}
	if(it_a != alpha_S[j].end() && it_b != beta_S[j].end()){
		const double prob = exp(it_a->second + it_b->second - logZ);
		cerr << "  pair_prob=" << prob << endl;
	}
	cerr << "  logZ=" << logZ << endl;

	if(i - 1 >= 0 && j + 1 < seq_n){
		const double loop_energy = _energy->energy_loop(i - 1, j + 1, i, j);
		cerr << "  energy_loop(i-1,j+1,i,j)=" << loop_energy << endl;
	}
}

void LinCapR::debug_prob(int idx) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_prob: idx out of range: " << idx << endl;
		return;
	}
	const double b = prob_B[idx];
	const double h = prob_H[idx];
	const double in = prob_I[idx];
	const double m = prob_M[idx];
	const double e = prob_E[idx];
	const double s = prob_S[idx];
	const double sum = b + h + in + m + e + s;
	cerr << "debug_prob i=" << idx
	     << " B=" << b
	     << " H=" << h
	     << " I=" << in
	     << " M=" << m
	     << " E=" << e
	     << " S=" << s
	     << " sum=" << sum << endl;
}

void LinCapR::debug_hairpin(int idx, int topn) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_hairpin: idx out of range: " << idx << endl;
		return;
	}
	if(topn <= 0) topn = 1;

	struct Item {
		int j;
		int k;
		double prob;
	};
	vector<Item> items;
	double sum_contrib = 0.0;
	const double logZ = alpha_O[seq_n - 1];

	for(int k = 0; k < seq_n; k++){
		for(const auto [j, score] : beta_SE[k]){
			if(idx < j || idx > k) continue;
			const int outer_i = j - 1;
			const int outer_j = k + 1;
			if(outer_i < 0 || outer_j >= seq_n) continue;
			const double new_score = exp(score - _energy->energy_hairpin(outer_i, outer_j) / _energy->kT() - logZ);
			items.push_back({j, k, new_score});
			sum_contrib += new_score;
		}
	}

	sort(items.begin(), items.end(), [](const Item& a, const Item& b){
		return a.prob > b.prob;
	});

	const auto base_at = [&](int pos) -> char {
		if(pos < 0 || pos >= seq_n) return 'N';
		return seq[pos];
	};

	cerr << "debug_hairpin i=" << idx << " candidates=" << items.size() << endl;
	cerr << "  sum_contrib=" << sum_contrib
	     << " prob_H[i]=" << prob_H[idx]
	     << " normalize_profiles=" << (normalize_profiles ? 1 : 0) << endl;
	const int limit = min(topn, static_cast<int>(items.size()));
	for(int t = 0; t < limit; t++){
		const auto& it = items[t];
		const int outer_i = it.j - 1;
		const int outer_j = it.k + 1;
		const int unp_l = outer_i + 1;
		const int unp_r = outer_j - 1;
		const bool outer_pairable = can_pair(outer_i, outer_j);
		const bool outer_in_alpha = lcr::dp::contains(alpha_S, outer_i, outer_j);
		const double loop_energy = _energy->energy_hairpin(outer_i, outer_j);
		const int raw_loop_len = outer_j - outer_i - 1;
		cerr << "  (j,k)=(" << it.j << "," << it.k << ")"
		     << " prob=" << it.prob
		     << " raw_outer=(" << outer_i << "," << outer_j << ")"
		     << " raw_unpaired=(" << unp_l << "," << unp_r << ")"
		     << " bases=" << base_at(outer_i) << "," << base_at(outer_j)
		     << " outer_can_pair=" << outer_pairable
		     << " outer_in_alpha=" << outer_in_alpha
		     << " raw_loop_len=" << raw_loop_len
		     << " hairpinE=" << loop_energy
		     << endl;
	}
}

void LinCapR::debug_external(int idx, int topn) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_external: idx out of range: " << idx << endl;
		return;
	}
	if(topn <= 0) topn = 1;
	const double logZ = alpha_O[seq_n - 1];
	const int left = idx - 1;
	const int right = idx + 1;
	const double alpha_left = (left >= 0 ? alpha_O[left] : 0.0);
	const double beta_right = (right < seq_n ? beta_O[right] : 0.0);
	const double kT = _energy->kT();
	const double unpaired_e_kcal = _energy->energy_external_unpaired(idx, idx);
	const double unpaired_e = unpaired_e_kcal / kT;
	double e_log = 0.0;
	if(idx == 0){
		e_log = beta_O[1] - unpaired_e - logZ;
	}else if(idx == seq_n - 1){
		e_log = alpha_O[seq_n - 2] - unpaired_e - logZ;
	}else{
		e_log = alpha_left + beta_right - unpaired_e - logZ;
	}
	cerr << "debug_external i=" << idx
	     << " alpha_O[i-1]=" << alpha_left
	     << " beta_O[i+1]=" << beta_right
	     << " kT=" << kT
	     << " unpairedE_kcal=" << unpaired_e_kcal
	     << " unpairedE=" << unpaired_e
	     << " logZ=" << logZ
	     << " E=" << exp(e_log)
	     << endl;

	struct Item {
		double logw;
		int i;
		int j;
		const char* label;
		double e_kcal;
		double e_scaled;
	};

	if(left >= 0){
		vector<Item> items;
		double total = -INF;
		if(left == 0){
			const double e_kcal = _energy->energy_external_unpaired(0, 0);
			const double e_scaled = e_kcal / kT;
			const double logw = -e_scaled;
			items.push_back({logw, 0, 0, "O_init", e_kcal, e_scaled});
			total = lcr::dp::logsumexp(total, logw);
		}else if(left - 1 >= 0){
			const double e_kcal = _energy->energy_external_unpaired(left, left);
			const double e_scaled = e_kcal / kT;
			const double logw = alpha_O[left - 1] - e_scaled;
			items.push_back({logw, left - 1, left, "O->O", e_kcal, e_scaled});
			total = lcr::dp::logsumexp(total, logw);
		}
		for(const auto [i, score] : alpha_S[left]){
			const double e_kcal = _energy->energy_external(i, left);
			const double e_scaled = e_kcal / kT;
			const double logw = (i - 1 >= 0 ? alpha_O[i - 1] : 0.0)
				+ score
				- e_scaled;
			items.push_back({logw, i, left, "O+S", e_kcal, e_scaled});
			total = lcr::dp::logsumexp(total, logw);
		}
		sort(items.begin(), items.end(), [](const Item& a, const Item& b){
			return a.logw > b.logw;
		});
		cerr << "  alpha_O[" << left << "] logsum=" << total
		     << " diff=" << (alpha_left - total)
		     << " items=" << items.size() << endl;
		const int limit = min(topn, static_cast<int>(items.size()));
		for(int t = 0; t < limit; t++){
			const auto& it = items[t];
			const double w = exp(it.logw - alpha_left);
			cerr << "    " << it.label
			     << " (i,j)=(" << it.i << "," << it.j << ")"
			     << " logw=" << it.logw
			     << " e_kcal=" << it.e_kcal
			     << " e_scaled=" << it.e_scaled
			     << " weight=" << w
			     << endl;
		}
	}

	if(right < seq_n){
		vector<Item> items;
		double total = -INF;
		if(right == seq_n - 1){
			const double e_kcal = _energy->energy_external_unpaired(seq_n - 1, seq_n - 1);
			const double e_scaled = e_kcal / kT;
			const double logw = -e_scaled;
			items.push_back({logw, seq_n - 1, seq_n - 1, "O_init", e_kcal, e_scaled});
			total = lcr::dp::logsumexp(total, logw);
		}else if(right + 1 < seq_n){
			const double e_kcal = _energy->energy_external_unpaired(right + 1, right + 1);
			const double e_scaled = e_kcal / kT;
			const double logw = beta_O[right + 1] - e_scaled;
			items.push_back({logw, right + 1, right, "O->O", e_kcal, e_scaled});
			total = lcr::dp::logsumexp(total, logw);
		}
		for(int j = right; j < seq_n; j++){
			auto it = alpha_S[j].find(right);
			if(it == alpha_S[j].end()) continue;
			const double e_kcal = _energy->energy_external(right, j);
			const double e_scaled = e_kcal / kT;
			const double logw = it->second
				+ (j + 1 < seq_n ? beta_O[j + 1] : 0.0)
				- e_scaled;
			items.push_back({logw, right, j, "O+S", e_kcal, e_scaled});
			total = lcr::dp::logsumexp(total, logw);
		}
		sort(items.begin(), items.end(), [](const Item& a, const Item& b){
			return a.logw > b.logw;
		});
		cerr << "  beta_O[" << right << "] logsum=" << total
		     << " diff=" << (beta_right - total)
		     << " items=" << items.size() << endl;
		const int limit = min(topn, static_cast<int>(items.size()));
		for(int t = 0; t < limit; t++){
			const auto& it = items[t];
			const double w = exp(it.logw - beta_right);
			cerr << "    " << it.label
			     << " (i,j)=(" << it.i << "," << it.j << ")"
			     << " logw=" << it.logw
			     << " e_kcal=" << it.e_kcal
			     << " e_scaled=" << it.e_scaled
			     << " weight=" << w
			     << endl;
		}
	}
}

void LinCapR::debug_raccess_outer_branch(int i, int j) const {
#ifdef LINEARRACCESS_WITH_RACCESS
	lcr::debug_raccess_outer_branch_terms(seq, i, j);
#else
	static_cast<void>(i);
	static_cast<void>(j);
	throw std::runtime_error("Raccess backend is unavailable in this build");
#endif
}

void LinCapR::debug_internal(int idx, int topn) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_internal: idx out of range: " << idx << endl;
		return;
	}
	if(topn <= 0) topn = 1;

	struct Item {
		int j;
		int k;
		int p;
		int q;
		double prob;
		const char* side;
	};
	vector<Item> items;
	const double logZ = alpha_O[seq_n - 1];

	for(int k = 0; k < seq_n; k++){
		for(const auto [j, score] : beta_SE[k]){
			for(int p = j; p <= min(j + MAXLOOP, k - 1); p++){
				for(int q = k; q >= p + TURN + 1 && (p - j) + (k - q) <= MAXLOOP; q--){
					if((p == j && q == k) || !lcr::dp::contains(alpha_S, p, q)) continue;
					auto it_a = alpha_S[q].find(p);
					if(it_a == alpha_S[q].end()) continue;
					const Float new_score = exp(score + it_a->second - _energy->energy_loop(j - 1, k + 1, p, q) / _energy->kT() - logZ);

					const bool left_internal = (q != k) && (idx >= j && idx <= (p - 1));
					const bool right_internal = (p != j) && (idx >= (q + 1) && idx <= k);
					if(left_internal){
						items.push_back({j, k, p, q, new_score, "left"});
					}
					if(right_internal){
						items.push_back({j, k, p, q, new_score, "right"});
					}
				}
			}
		}
	}

	sort(items.begin(), items.end(), [](const Item& a, const Item& b){
		return a.prob > b.prob;
	});

	const auto base_at = [&](int pos) -> char {
		if(pos < 0 || pos >= seq_n) return 'N';
		return seq[pos];
	};

	cerr << "debug_internal i=" << idx << " candidates=" << items.size() << endl;
	const int limit = min(topn, static_cast<int>(items.size()));
	for(int t = 0; t < limit; t++){
		const auto& it = items[t];
		const int outer_i = it.j - 1;
		const int outer_j = it.k + 1;
		const int r_outer_i = outer_i - 1;
		const int r_outer_j = outer_j;
		const int r_inner_i = it.p - 1;
		const int r_inner_j = it.q;
		const int left_len = it.p - it.j;
		const int right_len = it.k - it.q;
		const bool outer_pairable = (0 <= outer_i && outer_j < seq_n) ? can_pair(outer_i, outer_j) : false;
		const bool outer_in_alpha = (0 <= outer_i && outer_j < seq_n) ? lcr::dp::contains(alpha_S, outer_i, outer_j) : false;
		const double loop_energy = (outer_i >= 0 && outer_j < seq_n)
			? _energy->energy_loop(outer_i, outer_j, it.p, it.q)
			: 0.0;
		cerr << "  " << it.side
		     << " (j,k,p,q)=(" << it.j << "," << it.k << "," << it.p << "," << it.q << ")"
		     << " prob=" << it.prob
		     << " outer=(" << outer_i << "," << outer_j << ")"
		     << " inner=(" << it.p << "," << it.q << ")"
		     << " raccess_outer=(" << r_outer_i << "," << r_outer_j << ")"
		     << " raccess_inner=(" << r_inner_i << "," << r_inner_j << ")"
		     << " len=(" << left_len << "," << right_len << ")"
		     << " bases outer=" << base_at(outer_i) << "," << base_at(outer_j)
		     << " inner=" << base_at(it.p) << "," << base_at(it.q)
		     << " raccess_bases outer=" << base_at(r_outer_i) << "," << base_at(r_outer_j)
		     << " inner=" << base_at(r_inner_i) << "," << base_at(r_inner_j)
		     << " outer_can_pair=" << outer_pairable
		     << " outer_in_alpha=" << outer_in_alpha
		     << " loopE=" << loop_energy
		     << endl;
	}
}

void LinCapR::debug_multi_unpaired(int i, int j) const {
	if(i < 0 || j < 0 || i >= seq_n || j >= seq_n || i > j){
		cerr << "debug_multi_unpaired: invalid range: " << i << "," << j << endl;
		return;
	}
	const int len_closed = j - i + 1;
	const int len_half = j - i;
	const double energy = _energy->energy_multi_unpaired(i, j);
	cerr << "debug_multi_unpaired (i,j)=(" << i << "," << j << ")"
	     << " len_closed=" << len_closed
	     << " len_half=" << len_half
	     << " energy=" << energy
	     << endl;
}

void LinCapR::debug_multi_prob(int idx, int topn) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_multi_prob: idx out of range: " << idx << endl;
		return;
	}
	if(topn <= 0) topn = 1;
	const Float logZ = alpha_O[seq_n - 1];
	double from_mb = 0.0;
	double from_s = 0.0;

	struct Item {
		int j;
		int q;
		int k;
		double prob;
	};
	vector<Item> items;

	// Contribution from alpha_MB / beta_M (calc_profile M part 1)
	for(int k = 0; k < seq_n; k++){
		for(const auto [p, score] : alpha_MB[k]){
			const int j_min = max(0, p - c_multi);
			for(int j = p - 1; j >= j_min; j--){
				if(!lcr::dp::contains(beta_M, j, k)) continue;
				if(idx < j || idx > (p - 1)) continue;
				auto it_b = beta_M[k].find(j);
				if(it_b == beta_M[k].end()) continue;
				const Float new_score = exp(score + it_b->second
					- _energy->energy_multi_unpaired(j, p - 1) / _energy->kT()
					- logZ);
				from_mb += new_score;
			}
		}
	}

	// Contribution from alpha_S / beta_M2 (calc_profile M part 2)
	for(int q = 0; q < seq_n; q++){
		for(const auto [j, score] : alpha_S[q]){
			const int k_max = min(seq_n - 1, q + c_multi);
			for(int k = q + 1; k <= k_max; k++){
				if(!lcr::dp::contains(beta_M2, j, k)) continue;
				if(idx < (q + 1) || idx > k) continue;
				auto it_b = beta_M2[k].find(j);
				if(it_b == beta_M2[k].end()) continue;
				const Float new_score = exp(score + it_b->second
					- (_energy->energy_multi_bif(j, q) + _energy->energy_multi_unpaired(q + 1, k)) / _energy->kT()
					- logZ);
				from_s += new_score;
				items.push_back({j, q, k, new_score});
			}
		}
	}

	const double total = from_mb + from_s;
	sort(items.begin(), items.end(), [](const Item& a, const Item& b){
		return a.prob > b.prob;
	});

	cerr << "debug_multi_prob i=" << idx
	     << " M_total=" << prob_M[idx]
	     << " M_from_mb=" << from_mb
	     << " M_from_s=" << from_s
	     << " M_sum=" << total
	     << endl;

	const int limit = min(topn, static_cast<int>(items.size()));
	for(int t = 0; t < limit; t++){
		const auto& it = items[t];
		cerr << "  M2 (j,q,k)=(" << it.j << "," << it.q << "," << it.k << ")"
		     << " prob=" << it.prob
		     << " range=(" << (it.q + 1) << "," << it.k << ")"
		     << endl;
	}
}

void LinCapR::debug_dp_dump(const string& state, int i0, int j0, int i1, int j1) const {
	if(seq_n <= 0) return;
	string key = state;
	for(char& c : key) c = static_cast<char>(tolower(static_cast<unsigned char>(c)));
	const Table* alpha = nullptr;
	const Table* beta = nullptr;
	const char* label = nullptr;
	if(key == "s" || key == "stem"){
		alpha = &alpha_S;
		beta = &beta_S;
		label = "S";
	}else if(key == "se" || key == "stem_end" || key == "stemend"){
		alpha = &alpha_SE;
		beta = &beta_SE;
		label = "SE";
	}else{
		cerr << "debug_dp_dump: unknown state: " << state << endl;
		return;
	}
	int lo_i = std::min(i0, i1);
	int hi_i = std::max(i0, i1);
	int lo_j = std::min(j0, j1);
	int hi_j = std::max(j0, j1);
	lo_i = std::max(0, lo_i);
	hi_i = std::min(seq_n - 1, hi_i);
	lo_j = std::max(0, lo_j);
	hi_j = std::min(seq_n - 1, hi_j);
	for(int i = lo_i; i <= hi_i; i++){
		for(int j = std::max(i, lo_j); j <= hi_j; j++){
			const Table& a_tbl = *alpha;
			const Table& b_tbl = *beta;
			const auto it_a = a_tbl[j].find(i);
			const auto it_b = b_tbl[j].find(i);
			const double a = (it_a == a_tbl[j].end() ? -INF : it_a->second);
			const double b = (it_b == b_tbl[j].end() ? -INF : it_b->second);
			cerr << "debug_dp_dump state=" << label
			     << " raw_pair=(" << i << "," << j << ")"
			     << " alpha=" << a
			     << " beta=" << b
			     << endl;
		}
	}
}

void LinCapR::debug_outer_dp_dump(int i0, int i1) const {
	if(seq_n <= 0){
		cerr << "debug_outer_dp: seq_n<=0 (" << seq_n << ")" << endl;
		return;
	}
	int lo = std::min(i0, i1);
	int hi = std::max(i0, i1);
	lo = std::max(0, lo);
	hi = std::min(seq_n - 1, hi);
	const double logZ = alpha_O[seq_n - 1];
	for(int i = lo; i <= hi; i++){
		const double alpha = alpha_O[i];
		const double beta = beta_O[i];
		const double beta_norm = beta - logZ;
		cerr << "debug_outer_dp idx=" << i
		     << " alpha_O=" << alpha
		     << " beta_O=" << beta
		     << " beta_norm=" << beta_norm
		     << " logZ=" << logZ
		     << endl;
	}
}
