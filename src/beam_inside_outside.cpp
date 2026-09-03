#include "beam_inside_outside.hpp"
#include "beam_prune.hpp"
#include "linearraccess/dp_table_api.hpp"
#include "linearraccess/seq_utils.hpp"
#include "energy_linearcapr.hpp"
#ifdef LINEARRACCESS_WITH_RACCESS
#include "energy_raccess.hpp"
#endif

#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cctype>
#include <stdexcept>

LinCapR::LinCapR(int beam_size, energy::Model model, EnergyEngine engine,
                 bool raccess_use_lincapr_external,
                 bool normalize_profiles, Float normalize_warn_eps, int c_multi,
                 int c_hairpin)
	: params(energy::get_params(model)),
	  beam_size(beam_size),
	  c_multi(c_multi),
	  c_hairpin(c_hairpin),
	  normalize_profiles(normalize_profiles),
	  normalize_warn_eps(normalize_warn_eps) {
	if(beam_size < 0){
		throw std::invalid_argument("beam size must be non-negative");
	}
	if(c_multi < 0){
		throw std::invalid_argument("c_multi must be non-negative");
	}
	if(c_hairpin < 0){
		throw std::invalid_argument("c_hairpin must be non-negative");
	}
	if(normalize_warn_eps < 0 || !std::isfinite(normalize_warn_eps)){
		throw std::invalid_argument(
		    "normalization warning tolerance must be finite and non-negative");
	}
	if(raccess_use_lincapr_external && engine != EnergyEngine::Raccess){
		throw std::invalid_argument(
		    "raccess_use_lincapr_external requires the Raccess engine");
	}
	switch (engine) {
	case EnergyEngine::Raccess:
#ifdef LINEARRACCESS_WITH_RACCESS
		_energy = lcr::make_raccess_energy(&params, raccess_use_lincapr_external);
#else
		throw std::invalid_argument(
		    "Raccess backend is unavailable in this build; rebuild with "
		    "LINEARRACCESS_WITH_RACCESS=ON and an authorized Raccess source tree");
#endif
		break;
	case EnergyEngine::LinearCapR:
		_energy.reset(new lcr::LinearCapREnergyModel(params));
		break;
	default:
		throw std::invalid_argument("unknown energy engine");
	}
	if(params.use_fast_logsumexp) set_logsumexp_fast_mode();
	else set_logsumexp_legacy_mode();
}

LinCapR::~LinCapR() {
}

void LinCapR::set_debug_hairpin_build(bool enable) {
	debug_hairpin_build_log = enable;
}

void LinCapR::set_debug_se(int i, int j, int max_hits) {
	debug_se_log = true;
	debug_se_i = i;
	debug_se_j = j;
	debug_se_hits = 0;
	debug_se_max_hits = max_hits;
}

void LinCapR::set_debug_beam(int j0, int j1, int interval, bool stop_after) {
	debug_beam_log = true;
	debug_beam_j0 = j0;
	debug_beam_j1 = j1;
	debug_beam_interval = max(1, interval);
	debug_beam_stop = stop_after;
	debug_beam_early_stop = false;
}

void LinCapR::set_debug_multi_energy(int i, int j) {
	debug_multi_energy = true;
	debug_multi_i = i;
	debug_multi_j = j;
	debug_multi_reported = false;
}

void LinCapR::set_debug_outer_alpha(int max_hits) {
	debug_outer_alpha_log = true;
	debug_outer_alpha_hits = 0;
	debug_outer_alpha_max = max(1, max_hits);
}

void LinCapR::set_debug_outer_dp_inside(int i0, int i1) {
	debug_outer_dp_inside = true;
	debug_outer_dp_inside_i0 = i0;
	debug_outer_dp_inside_i1 = i1;
}

void LinCapR::release_energy(bool leak) {
	if(leak){
		_energy.release();
		cerr << "debug_release_energy: leaked energy model" << endl;
		return;
	}
	_energy.reset();
	cerr << "debug_release_energy: reset energy model" << endl;
}

void LinCapR::set_debug_skip_outside(bool enable) {
	debug_skip_outside = enable;
}

void LinCapR::set_debug_skip_profile(bool enable) {
	debug_skip_profile = enable;
}

double LinCapR::logW_all_unpaired() const {
	if(seq_n <= 0) return 0.0;
	double sum = 0.0;
	for(int i = 0; i < seq_n; i++){
		sum += _energy->energy_external_unpaired(i, i) / _energy->kT();
	}
	return -sum;
}


// prune top-k states
Float LinCapR::prune(Map<int, Float> &states) const{
	return lcr::beam::prune_states(states, beam_size,
				       [this](const int i, const Float score) {
					 return (i >= 1 ? alpha_O[i - 1] : Float(0)) + score;
				       });
}


// clear temp tables & profiles
void LinCapR::clear(){
	seq = "";
	seq_int.clear();
	seq_n = 0;
	alpha_O.clear();
	beta_O.clear();

	for(int i = 0; i < NTABLES; i++){
		alphas[i]->clear();
		betas[i]->clear();
	}

	for(int i = 0; i < NPROBS; i++) probs[i]->clear();
}


// returns free energy of ensemble in kcal/mol
Float LinCapR::get_energy_ensemble() const{
	const Float logZ = alpha_O[seq_n - 1];
	if (dynamic_cast<const lcr::LinearCapREnergyModel*>(_energy.get()) != nullptr) {
		return (logZ * -(params.temperature + params.k0) * params.gas_constant) / 1000;
	}
	return -_energy->kT() * logZ;
}

Float LinCapR::get_logZ() const{
	return alpha_O[seq_n - 1];
}

Float LinCapR::get_kT() const{
	return _energy->kT();
}

const vector<Float>& LinCapR::get_prob_stem() const{
	return prob_S;
}

const vector<Float>& LinCapR::get_prob_bulge() const{
	return prob_B;
}

const vector<Float>& LinCapR::get_prob_exterior() const{
	return prob_E;
}

const vector<Float>& LinCapR::get_prob_hairpin() const{
	return prob_H;
}

const vector<Float>& LinCapR::get_prob_internal() const{
	return prob_I;
}

const vector<Float>& LinCapR::get_prob_multiloop() const{
	return prob_M;
}


// calc structural profile
void LinCapR::run(const string &seq){
	initialize(lcr::seq::normalize_sequence(seq));
	if(debug_multi_energy && !debug_multi_reported){
		const int i = debug_multi_i;
		const int j = debug_multi_j;
		if(i >= 0 && j >= 0 && i <= j && j < seq_n){
			const double kT = _energy->kT();
			const double e_ext_unp = _energy->energy_external_unpaired(i, j);
			const double e_ext_br = (i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1))
				? _energy->energy_external(i - 1, j + 1)
				: NAN;
			const double e_bif = _energy->energy_multi_bif(i, j);
			const double e_close = (i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1))
				? _energy->energy_multi_closing(i - 1, j + 1)
				: NAN;
			const double e_unp = (i + 1 <= j) ? _energy->energy_multi_unpaired(i + 1, j) : 0.0;
			cerr << "debug_multi_energy: i=" << i << " j=" << j
			     << " kT=" << kT
			     << " ext_unpaired=" << e_ext_unp
			     << " ext_branch=" << e_ext_br
			     << " multi_bif=" << e_bif
			     << " multi_closing=" << e_close
			     << " multi_unpaired=" << e_unp
			     << endl;
			debug_multi_reported = true;
		}else{
			cerr << "debug_multi_energy: invalid i,j: i=" << i << " j=" << j
			     << " seq_n=" << seq_n << endl;
			debug_multi_reported = true;
		}
	}
	calc_inside();
	if(debug_outer_dp_inside){
		debug_outer_dp_dump(debug_outer_dp_inside_i0, debug_outer_dp_inside_i1);
	}
	if(debug_beam_early_stop){
		cerr << "debug_beam: early stop after calc_inside()" << endl;
		return;
	}
	if(debug_skip_outside){
		cerr << "debug_skip_outside: skip outside/profile" << endl;
		return;
	}
	calc_outside();
	if(debug_skip_profile){
		cerr << "debug_skip_profile: skip profile" << endl;
		return;
	}
	calc_profile();
}


// initialize
void LinCapR::initialize(const string &seq){
	this->seq = seq;

	// integerize sequence
	seq_n = seq.length();
	seq_int.resize(seq_n);
	for(int i = 0; i < seq_n; i++){
		seq_int[i] = lcr::seq::base_to_num(seq[i]);
	}
	_energy->set_sequence(seq, seq_int);

	// prepare DP tables
	alphas[0] = &alpha_S;
	alphas[1] = &alpha_SE;
	alphas[2] = &alpha_M;
	alphas[3] = &alpha_MB;
	alphas[4] = &alpha_M1;
	alphas[5] = &alpha_M2;

	betas[0] = &beta_S;
	betas[1] = &beta_SE;
	betas[2] = &beta_M;
	betas[3] = &beta_MB;
	betas[4] = &beta_M1;
	betas[5] = &beta_M2;

	alpha_O.assign(seq_n, -INF);
	beta_O.assign(seq_n, -INF);
	for(int i = 0; i < NTABLES; i++){
		alphas[i]->assign(seq_n, {});
		betas[i]->assign(seq_n, {});
		// google hash
		// for(int j = 0; j < seq_n; j++){
		// 	alphas[i]->at(j).set_empty_key(-1);
		// 	alphas[i]->at(j).set_deleted_key(-2);
		// 	betas[i]->at(j).set_empty_key(-1);
		// 	betas[i]->at(j).set_deleted_key(-2);
		// }
	}

	// prepare prob vectors
	probs[0] = &prob_B;
	probs[1] = &prob_E;
	probs[2] = &prob_H;
	probs[3] = &prob_I;
	probs[4] = &prob_M;
	probs[5] = &prob_S;

	for(int i = 0; i < NPROBS; i++) probs[i]->assign(seq_n, 0.0);

	// calc next pair index
	lcr::seq::build_next_pair(seq_int, next_pair);
}


// calc inside variables
void LinCapR::calc_inside(){
	alpha_O[0] = - _energy->energy_external_unpaired(0, 0) / _energy->kT();
	if(debug_outer_alpha_log){
		cerr << "debug_outer_alpha_state: phase=init"
		     << " alpha_O0=" << alpha_O[0]
		     << " alpha_O_last=" << (seq_n > 0 ? alpha_O[seq_n - 1] : 0.0)
		     << endl;
	}
	const auto base_at = [&](int idx) -> char {
		if(idx < 0 || idx >= seq_n) return 'N';
		return seq[idx];
	};
	const auto log_se = [&](const char* label, int i, int j, int outer_i, int outer_j,
	                        int inner_i, int inner_j, double base, double dsc, double logw) {
		if(!debug_se_log) return;
		if(debug_se_hits >= debug_se_max_hits) return;
		if(i != debug_se_i || j != debug_se_j) return;
		cerr << "debug_se: label=" << label
		     << " target=(" << i << "," << j << ")"
		     << " outer=(" << outer_i << "," << outer_j << ")"
		     << " inner=(" << inner_i << "," << inner_j << ")"
		     << " bases=" << base_at(outer_i) << "," << base_at(outer_j)
		     << " base=" << base
		     << " dsc=" << dsc
		     << " logw=" << logw
		     << endl;
		debug_se_hits++;
	};

	for(int j = 0; j < seq_n; j++){
		if(debug_multi_energy && !debug_multi_reported && j == debug_multi_j){
			const int i = debug_multi_i;
			if(i >= 0 && i < seq_n && i <= j){
				const double kT = _energy->kT();
				const double e_bif = _energy->energy_multi_bif(i, j);
				const double e_close = (i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1))
					? _energy->energy_multi_closing(i - 1, j + 1)
					: NAN;
				const double e_unp = (i + 1 <= j) ? _energy->energy_multi_unpaired(i + 1, j) : 0.0;
				cerr << "debug_multi_energy: i=" << i << " j=" << j
				     << " kT=" << kT
				     << " multi_bif=" << e_bif
				     << " multi_closing=" << e_close
				     << " multi_unpaired=" << e_unp
				     << endl;
				debug_multi_reported = true;
			}
		}
		const bool log_beam = debug_beam_log && j >= debug_beam_j0 && j <= debug_beam_j1
		                      && ((j - debug_beam_j0) % debug_beam_interval == 0);
		int upd_m2 = 0;
		int upd_mb = 0;
		int upd_m1 = 0;
		int upd_m = 0;
		auto table_stats = [](const Map<int, Float>& table) -> tuple<int, int, int, double, double> {
			if(table.empty()) return {0, 0, 0, 0.0, 0.0};
			int min_i = std::numeric_limits<int>::max();
			int max_i = std::numeric_limits<int>::min();
			double min_s = std::numeric_limits<double>::infinity();
			double max_s = -std::numeric_limits<double>::infinity();
			for(const auto [i, score] : table){
				if(i < min_i) min_i = i;
				if(i > max_i) max_i = i;
				if(score < min_s) min_s = score;
				if(score > max_s) max_s = score;
			}
			return {static_cast<int>(table.size()), min_i, max_i, min_s, max_s};
		};
		auto update_m2 = [&](int i, int jj, Float score) {
			++upd_m2;
			lcr::dp::update_sum(alpha_M2, i, jj, score);
		};
		auto update_mb = [&](int i, int jj, Float score) {
			++upd_mb;
			lcr::dp::update_sum(alpha_MB, i, jj, score);
		};
		auto update_m1 = [&](int i, int jj, Float score) {
			++upd_m1;
			lcr::dp::update_sum(alpha_M1, i, jj, score);
		};
		auto update_m = [&](int i, int jj, Float score) {
			++upd_m;
			lcr::dp::update_sum(alpha_M, i, jj, score);
		};
		auto log_prune = [&](const char* label, Map<int, Float>& table) {
			if(!log_beam) return;
			const int before = static_cast<int>(table.size());
			const Float thresh = prune(table);
			const int after = static_cast<int>(table.size());
			cerr << "debug_beam: j=" << j
			     << " table=" << label
			     << " before=" << before
			     << " after=" << after
			     << " threshold=" << thresh
			     << endl;
		};
		// S
		if(log_beam) log_prune("S", alpha_S[j]);
		else prune(alpha_S[j]);
		for(const auto [i, score] : alpha_S[j]){
			// S -> S
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				lcr::dp::update_sum(alpha_S, i - 1, j + 1, score - _energy->energy_loop(i - 1, j + 1, i, j) / _energy->kT());
			}

			// M2 -> S
			for(int n = 0; n <= c_multi && j + n < seq_n; n++){
				update_m2(i, j + n, score - (_energy->energy_multi_bif(i, j) + _energy->energy_multi_unpaired(j + 1, j + n)) / _energy->kT());
			}

			// SE -> S: p..i..j..q, [p - 1, q] can be pair
			for(int p = i; i - p <= MAXLOOP && p >= 1; p--){
				for(int q = next_pair[seq_int[p - 1]][j + 1]; q < seq_n && (q - j - 1) + (i - p) <= MAXLOOP; q = next_pair[seq_int[p - 1]][q + 1]){
					if((p == i && q == j + 1)) continue;
					const double dsc = _energy->energy_loop(p - 1, q, i, j) / _energy->kT();
					const double logw = score - dsc;
					log_se("SE_from_S_loop", p, q - 1, p - 1, q, i, j, score, dsc, logw);
					lcr::dp::update_sum(alpha_SE, p, q - 1, logw);
				}
			}

			// O -> O + S
			const double base = (i - 1 >= 0 ? alpha_O[i - 1] : 0.0);
			const double dsc = _energy->energy_external(i, j) / _energy->kT();
			const double logw = base + score - dsc;
			if(debug_outer_alpha_log && debug_outer_alpha_hits < debug_outer_alpha_max){
				cerr << "debug_outer_alpha: i=" << i
				     << " j=" << j
				     << " base=" << base
				     << " score=" << score
				     << " dsc=" << dsc
				     << " logw=" << logw
				     << endl;
				debug_outer_alpha_hits++;
			}
			lcr::dp::update_sum(alpha_O, j, logw);
		}

		// M2
		if(log_beam) log_prune("M2", alpha_M2[j]);
		else prune(alpha_M2[j]);
		for(const auto [i, score] : alpha_M2[j]){
			// M1 -> M2
			update_m1(i, j, score);

			// MB -> M1 + M2
			if(i - 1 >= 0){
				for(const auto [k, score_m1] : alpha_M1[i - 1]){
					update_mb(k, j, score_m1 + score);
				}
			}
		}

		// MB
		if(log_beam) log_prune("MB", alpha_MB[j]);
		else prune(alpha_MB[j]);
		for(const auto [i, score] : alpha_MB[j]){
			// M1 -> MB
			update_m1(i, j, score);

			// M -> MB
			for(int n = 0; n <= c_multi && i - n >= 0; n++){
				update_m(i - n, j, score);
			}
		}

		// M1
		if(log_beam) log_prune("M1", alpha_M1[j]);
		else prune(alpha_M1[j]);

		// M
		if(log_beam) log_prune("M", alpha_M[j]);
		else prune(alpha_M[j]);
		for(const auto [i, score] : alpha_M[j]){
			// SE -> M
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				const double dsc = _energy->energy_multi_closing(i - 1, j + 1) / _energy->kT();
				const double logw = score - dsc;
				log_se("SE_from_M_close", i, j, i - 1, j + 1, i, j, score, dsc, logw);
				lcr::dp::update_sum(alpha_SE, i, j, logw);
			}
		}

		// SE -> (Hairpin)
		for(int n = TURN; n <= std::min(c_hairpin, j + 1); n++){
			const int i = j - n + 1;
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				if(debug_hairpin_build_log){
					cerr << "debug_hairpin_build raw_outer=(" << (i - 1) << "," << (j + 1)
					     << ") raw_loop_len=" << (j - i)
					     << " inner=(" << i << "," << j << ")"
					     << " bases=" << seq[i - 1] << "," << seq[j + 1]
					     << endl;
				}
				const double dsc = _energy->energy_hairpin(i - 1, j + 1) / _energy->kT();
				const double logw = -dsc;
				log_se("SE_from_hairpin", i, j, i - 1, j + 1, i, j, 0.0, dsc, logw);
				lcr::dp::update_sum(alpha_SE, i, j, logw);
			}
		}

		// SE
		if(log_beam) log_prune("SE", alpha_SE[j]);
		else prune(alpha_SE[j]);
		for(const auto [i, score] : alpha_SE[j]){
			// S -> SE
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				lcr::dp::update_sum(alpha_S, i - 1, j + 1, score);
			}
		}
		if(log_beam){
			cerr << "debug_beam_updates: j=" << j
			     << " M2=" << upd_m2
			     << " MB=" << upd_mb
			     << " M1=" << upd_m1
			     << " M=" << upd_m
			     << endl;
			const auto [sz_m2, min_m2, max_m2, min_s_m2, max_s_m2] = table_stats(alpha_M2[j]);
			const auto [sz_mb, min_mb, max_mb, min_s_mb, max_s_mb] = table_stats(alpha_MB[j]);
			const auto [sz_m1, min_m1, max_m1, min_s_m1, max_s_m1] = table_stats(alpha_M1[j]);
			const auto [sz_m, min_m, max_m, min_s_m, max_s_m] = table_stats(alpha_M[j]);
			cerr << "debug_beam_stats: j=" << j
			     << " M2(size=" << sz_m2 << ",i=[" << min_m2 << "," << max_m2 << "],score=[" << min_s_m2 << "," << max_s_m2 << "])"
			     << " MB(size=" << sz_mb << ",i=[" << min_mb << "," << max_mb << "],score=[" << min_s_mb << "," << max_s_mb << "])"
			     << " M1(size=" << sz_m1 << ",i=[" << min_m1 << "," << max_m1 << "],score=[" << min_s_m1 << "," << max_s_m1 << "])"
			     << " M(size=" << sz_m << ",i=[" << min_m << "," << max_m << "],score=[" << min_s_m << "," << max_s_m << "])"
			     << endl;
		}
		if(debug_beam_log && debug_beam_stop && j >= debug_beam_j1){
			debug_beam_early_stop = true;
			break;
		}

		// O -> O
		if(j + 1 < seq_n){
			lcr::dp::update_sum(alpha_O, j + 1, alpha_O[j] - _energy->energy_external_unpaired(j + 1, j + 1) / _energy->kT());
		}
	}
	if(debug_outer_alpha_log){
		int finite = 0;
		for(int i = 0; i < seq_n; i++){
			if(alpha_O[i] > -INF / 2) finite++;
		}
		cerr << "debug_outer_alpha_state: phase=after_inside"
		     << " alpha_O0=" << (seq_n > 0 ? alpha_O[0] : 0.0)
		     << " alpha_O_last=" << (seq_n > 0 ? alpha_O[seq_n - 1] : 0.0)
		     << " finite=" << finite
		     << endl;
	}
}


// calc outside variables
void LinCapR::calc_outside(){
	if(debug_outer_alpha_log){
		int finite = 0;
		for(int i = 0; i < seq_n; i++){
			if(alpha_O[i] > -INF / 2) finite++;
		}
		cerr << "debug_outer_alpha_state: phase=before_outside"
		     << " alpha_O0=" << (seq_n > 0 ? alpha_O[0] : 0.0)
		     << " alpha_O_last=" << (seq_n > 0 ? alpha_O[seq_n - 1] : 0.0)
		     << " finite=" << finite
		     << endl;
	}
	beta_O[seq_n - 1] = - _energy->energy_external_unpaired(seq_n - 1, seq_n - 1) / _energy->kT();
	for(int j = seq_n - 1; j >= 0; j--){
		// O
		// O -> O
		if  (j + 1 < seq_n){
			lcr::dp::update_sum(beta_O, j, (j + 1 < seq_n ? beta_O[j + 1] : 0) - _energy->energy_external_unpaired(j + 1, j + 1) / _energy->kT());
		}

		// O -> O + S
		for(const auto [i, score] : alpha_S[j]){
			lcr::dp::update_sum(beta_O, i, score + (j + 1 < seq_n ? beta_O[j + 1] : 0) - _energy->energy_external(i, j) / _energy->kT());
		}

		// SE
		for(const auto [i, _] : alpha_SE[j]){
			// S -> SE
			if(i - 1 >= 0 && j + 1 < seq_n){
				lcr::dp::update_sum(beta_SE, i, j, lcr::dp::get_value(beta_S, i - 1, j + 1));
			}
		}

		// M
		for(const auto [i, _] : alpha_M[j]){
			// SE -> M
			if(i - 1 >= 0 && j + 1 < seq_n){
				lcr::dp::update_sum(beta_M, i, j, lcr::dp::get_value(beta_SE, i, j) - _energy->energy_multi_closing(i - 1, j + 1) / _energy->kT());
			}
		}

		// MB
		for(const auto [i, _] : alpha_MB[j]){
			// M1 -> MB
			lcr::dp::update_sum(beta_MB, i, j, lcr::dp::get_value(beta_M1, i, j));

			// M -> MB
			for(int n = 0; n <= c_multi && i - n >= 0; n++){
				lcr::dp::update_sum(beta_MB, i, j, lcr::dp::get_value(beta_M, i - n, j));
			}
		}

		// M1, M2
		for(const auto [i, score_M2] : alpha_M2[j]){
			// M1 -> M2
			lcr::dp::update_sum(beta_M2, i, j, lcr::dp::get_value(beta_M1, i, j));

			// MB -> M1 + M2
			if(i - 1 < 0) continue;
			for(const auto [k, score_M1] : alpha_M1[i - 1]){
				lcr::dp::update_sum(beta_M1, k, i - 1, lcr::dp::get_value(beta_MB, k, j) + score_M2);
				lcr::dp::update_sum(beta_M2, i, j, lcr::dp::get_value(beta_MB, k, j) + score_M1);
			}
		}

		// S
		for(const auto [i, _] : alpha_S[j]){
			// O -> O + S
			lcr::dp::update_sum(beta_S, i, j, (i - 1 >= 0 ? alpha_O[i - 1] : 0) + (j + 1 < seq_n ? beta_O[j + 1] : 0) - _energy->energy_external(i, j) / _energy->kT());

			// SE -> S
			for(int p = i; i - p <= MAXLOOP && p >= 1; p--){
				for(int q = next_pair[seq_int[p - 1]][j + 1]; q < seq_n && (q - j - 1) + (i - p) <= MAXLOOP; q = next_pair[seq_int[p - 1]][q + 1]){
					if((p == i && q == j + 1)) continue;
					lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_SE, p, q - 1) - _energy->energy_loop(p - 1, q, i, j) / _energy->kT());
				}
			}

			// S -> S
			if(i - 1 >= 0 && j + 1 < seq_n){
				lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_S, i - 1, j + 1) - _energy->energy_loop(i - 1, j + 1, i, j) / _energy->kT());
			}

			// M2 -> S
			for(int n = 0; n <= c_multi && j + n < seq_n; n++){
				lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_M2, i, j + n) - (_energy->energy_multi_bif(i, j) + _energy->energy_multi_unpaired(j + 1, j + n)) / _energy->kT());
			}
		}
	}
	if(debug_outer_alpha_log){
		int finite = 0;
		for(int i = 0; i < seq_n; i++){
			if(alpha_O[i] > -INF / 2) finite++;
		}
		cerr << "debug_outer_alpha_state: phase=after_outside"
		     << " alpha_O0=" << (seq_n > 0 ? alpha_O[0] : 0.0)
		     << " alpha_O_last=" << (seq_n > 0 ? alpha_O[seq_n - 1] : 0.0)
		     << " finite=" << finite
		     << endl;
	}
}


// calc energy of hairpin loop [i, j]
