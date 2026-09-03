#include "energy_raccess.hpp"
#include "energy_linearcapr.hpp"

// Avoid macro conflicts (e.g. MAXLOOP) with Raccess headers.
#ifdef MAXLOOP
#pragma push_macro("MAXLOOP")
#undef MAXLOOP
#define LCR_RESTORE_MAXLOOP 1
#endif
#include "raccess/energy_model_api.hpp"
#include "raccess/prob_model.hpp"
#include "util/util.hpp"
#include <algorithm>
#include <cstdio>
#include <limits>
#ifdef LCR_RESTORE_MAXLOOP
#pragma pop_macro("MAXLOOP")
#undef LCR_RESTORE_MAXLOOP
#endif

namespace lcr {

class RaccessEnergyModelImpl : public EnergyApi {
public:
	static constexpr double kLengthFactor = 0.0;
	static constexpr bool kDebugExternalUnpairedZero = false;
	static constexpr bool kDebugExternalBranchZero = false;
	RaccessEnergyModelImpl(const energy::Params* params, bool use_linearcapr_external)
		: _sm(), _api(_sm), _use_linearcapr_external(use_linearcapr_external) {
		_api.initialize();
		// Choose the zero gauge.  Upstream Raccess uses -0.541728723 by
		// default, but applying the same per-nucleotide factor to paired and
		// unpaired exterior transitions leaves exact normalized probabilities
		// unchanged (and avoids biasing beam pruning scores here).
		_sm.set_param("external_unpaired", kLengthFactor);
		_sm.set_param("external_paired_length_factor", kLengthFactor);
		if(_use_linearcapr_external){
			if(params == nullptr){
				std::cerr << "RaccessEnergyModelImpl: linearcapr external requested without params; disabling override."
				          << std::endl;
				_use_linearcapr_external = false;
			}else{
				_external_linearcapr.reset(new LinearCapREnergyModel(*params));
			}
		}
	}
	~RaccessEnergyModelImpl() override {
	}
	void set_sequence(const std::string& seq, const std::vector<int>& seq_int) override {
		SM::Seq codes;
		codes.resize(seq.size());
		::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
		_api.set_seq(codes);
		if(_external_linearcapr){
			_external_linearcapr->set_sequence(seq, seq_int);
		}
	}
	double kT() const override { return Raccess::ScoreModelEnergy::RT_KCAL_MOL(); }
	Float energy_hairpin(int i, int j) const override {
		const int a = i + 1;
		const int b = j + 1;
		const double energy = _api.score_to_energy(_api.log_boltz_hairpin_closed(a, b));
		return energy * debug_kHairpinScale;
	}
	Float energy_loop(int i, int j, int p, int q) const override {
		// const int a = i + 1;
		// const int b = j + 1;
		// const int c = p + 1;
		// const int d = q + 1;
		// return _api.score_to_energy(_api.log_boltz_loop_closed(a, b, c, d));
		// Use DP-coordinate mapping to align with Raccess interior loop indexing.
		if (p == i + 1 && q == j - 1) {
			const int a = i + 1;
			const int b = j + 1;
			return _api.score_to_energy(_api.log_boltz_stack_closed(a, b));
		}
		const int dp_i = i + 1;
		const int dp_j = j;
		const int dp_p = p;
		const int dp_q = q + 1;
		return _api.score_to_energy(_api.log_boltz_interior(dp_i, dp_j, dp_p, dp_q));
	}
	Float energy_external(int i, int j) const override {
		if(_use_linearcapr_external && _external_linearcapr){
			return _external_linearcapr->energy_external(i, j);
		}
		if(kDebugExternalBranchZero){
			return 0.0;
		}
		const int a = i + 1;
		const int b = j + 1;
		return _api.score_to_energy(_api.log_boltz_outer_branch_closed(a, b));
	}
	Float energy_external_unpaired(int i, int j) const override {
		if(_use_linearcapr_external && _external_linearcapr){
			return _external_linearcapr->energy_external_unpaired(i, j);
		}
		if(kDebugExternalUnpairedZero){
			return 0.0;
		}
		const int len = (j - i + 1);
		return _api.score_to_energy(_api.log_boltz_outer_extend(0, len));
	}
	Float energy_multi_unpaired(int i, int j) const override {
		const int len = (j - i + 1);
		return _api.score_to_energy(_api.log_boltz_multi_extend(0, len));
	}
	Float energy_multi_closing(int i, int j) const override {
		const int a = i + 1;
		const int b = j + 1;
		return _api.score_to_energy(_api.log_boltz_multi_close_closed(a, b));
	}
	Float energy_multi_bif(int i, int j) const override {
		const int a = i + 1;
		const int b = j + 1;
		return _api.score_to_energy(_api.log_boltz_multi_open_closed(a, b));
	}

private:
	static constexpr double debug_kHairpinScale = 1.0;
	typedef Raccess::ScoreModelEnergy SM;
	SM _sm;
	Raccess::EnergyModelApi _api;
	bool _use_linearcapr_external = false;
	std::unique_ptr<LinearCapREnergyModel> _external_linearcapr;
};

std::unique_ptr<EnergyApi> make_raccess_energy(const energy::Params* params,
                                               bool use_linearcapr_external) {
	return std::unique_ptr<EnergyApi>(new RaccessEnergyModelImpl(params, use_linearcapr_external));
}

double raccess_rt_kcal_mol() {
	return Raccess::ScoreModelEnergy::RT_KCAL_MOL();
}

double raccess_length_factor() {
	return RaccessEnergyModelImpl::kLengthFactor;
}

std::vector<double> compute_raccess_unpaired_1(const std::string& seq, int max_span) {
	const auto result = compute_raccess_accessibility(seq, {1}, max_span, 0.0);
	return result.empty() ? std::vector<double>() : result.front();
}

std::vector<std::vector<double>> compute_raccess_accessibility(
		const std::string& seq, const std::vector<int>& access_lens, int max_span,
		double length_factor) {
	using SM = Raccess::ScoreModelEnergy;
	using PM = Raccess::ProbModel<SM>;
	SM sm;
	sm.initialize();
	sm.set_param("external_unpaired", length_factor);
	sm.set_param("external_paired_length_factor", length_factor);

	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	sm.set_seq(codes);

	PM pm;
	pm.set_score_model(sm);
	pm.set_max_span(max_span);
	pm.set_prob_thr(0);
	PM::VI acc_lens;
	for(const int length : access_lens){
		if(length > 0 && length <= static_cast<int>(seq.size())){
			acc_lens.push_back(length);
		}
	}
	if(acc_lens.empty()) return {};
	pm.set_acc_lens(acc_lens);

	std::vector<std::vector<double>> probabilities;
	probabilities.reserve(acc_lens.size());
	for(const int length : acc_lens){
		probabilities.emplace_back(seq.size() - length + 1,
		                           std::numeric_limits<double>::quiet_NaN());
	}
	auto block = [&](int i, int w, double energy) {
		auto it = std::find(acc_lens.begin(), acc_lens.end(), w);
		if(it == acc_lens.end()) return;
		const double logp = sm.energy_to_score(energy);
		const size_t length_idx = static_cast<size_t>(it - acc_lens.begin());
		probabilities[length_idx][i] = exp(logp);
	};
	pm.compute_prob(block);
	return probabilities;
}

void debug_raccess_multi_unpaired(const std::string& seq, int i, int j) {
	using SM = Raccess::ScoreModelEnergy;
	using Api = Raccess::EnergyModelApi;
	SM sm;
	sm.initialize();

	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	sm.set_seq(codes);

	Api api(sm);

	const int len_closed = j - i + 1;
	const int len_half = j - i;
	const double e_closed = api.score_to_energy(api.log_boltz_multi_extend(0, len_closed));
	const double e_half = api.score_to_energy(api.log_boltz_multi_extend(0, len_half));
	std::fprintf(stderr,
	             "raccess_multi_unpaired (i,j)=(%d,%d) len_closed=%d len_half=%d energy_closed=%g energy_half=%g\n",
	             i, j, len_closed, len_half, e_closed, e_half);
}

void debug_raccess_local(const std::string& seq, int i, int j, bool has_loop, int p, int q) {
	using SM = Raccess::ScoreModelEnergy;
	using Api = Raccess::EnergyModelApi;
	SM sm;
	sm.initialize();

	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	sm.set_seq(codes);

	Api api(sm);

	const int a = i + 1;
	const int b = j + 1;
	const int n = static_cast<int>(seq.size());

	const auto base_at = [&](int one_origin) -> char {
		if(one_origin < 1 || one_origin > n) return 'N';
		return seq[one_origin - 1];
	};

	const auto report = [&](const char* label, double via_closed, double direct) {
		const double diff = via_closed - direct;
		std::fprintf(stderr, "%s (%d,%d): closed=%g direct=%g diff=%g\n", label, i, j, via_closed, direct, diff);
	};

	report("hairpin",
	       api.score_to_energy(api.log_boltz_hairpin_closed(a, b)),
	       api.score_to_energy(api.log_boltz_hairpin(a, b - 1)));

	report("multi_close",
	       api.score_to_energy(api.log_boltz_multi_close_closed(a, b)),
	       api.score_to_energy(api.log_boltz_multi_close(a + 1, b - 1)));

	report("multi_open",
	       api.score_to_energy(api.log_boltz_multi_open_closed(a, b)),
	       api.score_to_energy(api.log_boltz_multi_open(a - 1, b)));

	report("external",
	       api.score_to_energy(api.log_boltz_outer_branch_closed(a, b)),
	       api.score_to_energy(api.log_boltz_outer_branch(a - 1, b)));

	if (has_loop) {
		const int c = p + 1;
		const int d = q + 1;
		const int dp_i = a + 1;
		const int dp_j = b - 1;
		std::fprintf(stderr, "loop mapping closed(a,b,c,d)=(%d,%d,%d,%d) dp(i,j,ip,jp)=(%d,%d,%d,%d)\n",
		             a, b, c, d, dp_i, dp_j, c, d);
		std::fprintf(stderr, "  closing pair (a,b) bases=%c,%c inner pair (c,d) bases=%c,%c\n",
		             base_at(a), base_at(b), base_at(c), base_at(d));
		std::fprintf(stderr,
		             "  padded bases: seq(i)=%c seq(i+1)=%c seq(j)=%c seq(j+1)=%c seq(ip)=%c seq(ip+1)=%c seq(jp)=%c seq(jp+1)=%c\n",
		             base_at(dp_i),
		             base_at(dp_i + 1),
		             base_at(dp_j),
		             base_at(dp_j + 1),
		             base_at(c),
		             base_at(c + 1),
		             base_at(d),
		             base_at(d + 1));
		const int li = c - a;
		const int lj = b - d;
		const double interior_score = sm.score_interior(dp_i, dp_j, c, d);
		const double interior_nuc_score = sm.score_interior_nuc(dp_i, dp_j, c, d);
		std::fprintf(stderr, "  interior lens (li,lj)=(%d,%d) score=%g nuc_score=%g\n",
		             li, lj, interior_score, interior_nuc_score);
		const double via_closed = api.score_to_energy(api.log_boltz_loop_closed(a, b, c, d));
		double direct = 0.0;
		if ((c == (a + 1)) && (d == (b - 1))) {
			direct = api.score_to_energy(api.log_boltz_stack(a, b + 1));
		} else {
			direct = api.score_to_energy(api.log_boltz_interior(a, b - 1, c - 1, d));
		}
		std::fprintf(stderr, "loop (%d,%d,%d,%d): closed=%g direct=%g diff=%g\n",
		             i, j, p, q, via_closed, direct, via_closed - direct);

		if ((c == (a + 1)) && (d == (b - 1))) {
			const double stack_ab1 = api.score_to_energy(api.log_boltz_stack(a, b + 1));
			const double stack_a1b = api.score_to_energy(api.log_boltz_stack(a - 1, b));
			std::fprintf(stderr, "stack variants (%d,%d): a,b+1=%g a-1,b=%g diff=%g\n",
			             i, j, stack_ab1, stack_a1b, stack_ab1 - stack_a1b);
		}
	}
}

void debug_raccess_outer_branch_terms(const std::string& seq, int i, int j) {
	using SM = Raccess::ScoreModelEnergy;
	using Api = Raccess::EnergyModelApi;
	if(i < 0 || j < 0 || i >= j || j >= static_cast<int>(seq.size())){
		std::fprintf(stderr, "debug_raccess_outer_branch_terms: invalid i,j (%d,%d) len=%zu\n",
		             i, j, seq.size());
		return;
	}
	SM sm;
	sm.initialize();

	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	sm.set_seq(codes);

	Api api(sm);
	const int a = i + 1;
	const int b = j + 1;

	const double branch_score = sm.score_outer_branch_nuc(a, b);
	const int idx_paired = sm.param_index("external_paired");
	const int idx_len = sm.param_index("external_paired_length_factor");
	const auto pv = sm.param_vector();
	const double external_paired = (idx_paired >= 0 ? pv[idx_paired] : 0.0);
	const double external_len = (idx_len >= 0 ? pv[idx_len] : 0.0);
	const double len_term = external_len * (j - i);
	const double total_score = branch_score + external_paired + len_term;
	const double api_score = api.log_boltz_outer_branch_closed(a, b);

	const auto to_energy = [&](double score) {
		return api.score_to_energy(score);
	};

	std::fprintf(stderr,
	             "debug_raccess_outer_branch_terms i=%d j=%d bases=%c,%c\n"
	             "  score_nuc=%g score_external_paired=%g score_external_len=%g score_len_term=%g\n"
	             "  score_total=%g score_api=%g score_diff=%g\n"
	             "  energy_nuc=%g energy_external_paired=%g energy_external_len=%g energy_len_term=%g\n"
	             "  energy_total=%g energy_api=%g energy_diff=%g\n",
	             i, j, seq[i], seq[j],
	             branch_score, external_paired, external_len, len_term,
	             total_score, api_score, total_score - api_score,
	             to_energy(branch_score), to_energy(external_paired), to_energy(external_len), to_energy(len_term),
	             to_energy(total_score), to_energy(api_score), to_energy(total_score) - to_energy(api_score));
}

double compute_raccess_logz(const std::string& seq, int max_span) {
	using SM = Raccess::ScoreModelEnergy;
	using PM = Raccess::ProbModel<SM>;
	SM sm;
	sm.initialize();

	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	sm.set_seq(codes);

	PM pm;
	pm.set_score_model(sm);
	pm.set_max_span(max_span);
	pm.set_prob_thr(0);
	PM::VI acc_lens;
	acc_lens.push_back(1);
	pm.set_acc_lens(acc_lens);

	auto no_op = [](int, int, double) {};
	pm.compute_prob(no_op);
	return pm.partition_coeff();
}

} // namespace lcr
