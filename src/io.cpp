#include "beam_inside_outside.hpp"
#include "io.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>

namespace lcr {
namespace io {

static void write_profile(std::ostream& os, const char* label, const std::vector<Float>& data) {
	os << label << " ";
	for(const auto& v : data) os << v << " ";
	os << std::endl;
}

static double energy_of_prob(double prob, double rt) {
	return (prob > 0.0)
		? -rt * std::log(prob)
		: std::numeric_limits<double>::infinity();
}

void loop_profile_io(std::ostream& os,
                     const std::string& seq_name,
                     const std::vector<Float>& prob_B,
                     const std::vector<Float>& prob_E,
                     const std::vector<Float>& prob_H,
                     const std::vector<Float>& prob_I,
                     const std::vector<Float>& prob_M,
                     const std::vector<Float>& prob_S) {
	os << ">" << seq_name << std::endl;
	write_profile(os, "Bulge", prob_B);
	write_profile(os, "Exterior", prob_E);
	write_profile(os, "Hairpin", prob_H);
	write_profile(os, "Internal", prob_I);
	write_profile(os, "Multiloop", prob_M);
	write_profile(os, "Stem", prob_S);
	os << std::endl;
}

void accessibility_raccess_io(std::ostream& os,
                              const std::string& seq_name,
                              int seq_len,
                              const std::vector<int>& access_lens,
                              const std::vector<std::vector<double>>& probs,
                              double rt) {
	const std::streamsize old_precision = os.precision();
	os << std::setprecision(std::numeric_limits<double>::max_digits10);
	os << ">" << seq_name << std::endl;
	if(access_lens.size() != probs.size()){
		os.precision(old_precision);
		return;
	}

	for(int pos = 0; pos < seq_len; pos++){
		bool wrote_any = false;
		os << pos;
		for(size_t idx = 0; idx < access_lens.size(); idx++){
			const int len = access_lens[idx];
			if(len <= 0) continue;
			const int last_start = seq_len - len;
			if(pos > last_start) continue;
			const auto& vec = probs[idx];
			if(pos >= static_cast<int>(vec.size())) continue;
			const double prob = vec[pos];
			const double energy = energy_of_prob(prob, rt);
			os << (wrote_any ? "" : "\t") << len << "," << energy << ";";
			wrote_any = true;
		}
		if(wrote_any){
			os << std::endl;
		}else{
			os << std::endl;
		}
	}
	os << std::endl;
	os.precision(old_precision);
}

void accessibility_probability_io(
        std::ostream& os,
        const std::string& seq_name,
        int seq_len,
        const std::vector<int>& access_lens,
        const std::vector<std::vector<double>>& probs) {
	const std::streamsize old_precision = os.precision();
	os << std::setprecision(std::numeric_limits<double>::max_digits10);
	os << ">" << seq_name << std::endl;
	os << "# value=probability" << std::endl;
	if(access_lens.size() != probs.size()){
		os.precision(old_precision);
		return;
	}

	for(int pos = 0; pos < seq_len; pos++){
		os << pos;
		bool wrote_any = false;
		for(size_t idx = 0; idx < access_lens.size(); idx++){
			const int len = access_lens[idx];
			if(len <= 0 || pos > seq_len - len) continue;
			const auto& vec = probs[idx];
			if(pos >= static_cast<int>(vec.size())) continue;
			os << (wrote_any ? "" : "\t") << len << "," << vec[pos] << ";";
			wrote_any = true;
		}
		os << std::endl;
	}
	os << std::endl;
	os.precision(old_precision);
}

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
                                     double rt) {
	const std::streamsize old_precision = os.precision();
	os << std::setprecision(std::numeric_limits<double>::max_digits10);
	os << ">" << seq_name << std::endl;
	const size_t n = access_lens.size();
	if(n == 0 || total_probs.size() != n || exterior_probs.size() != n ||
	   hairpin_probs.size() != n || bulge_probs.size() != n ||
	   internal_probs.size() != n || multiloop_probs.size() != n) {
		os << std::endl;
		os.precision(old_precision);
		return;
	}

	for(int pos = 0; pos < seq_len; pos++){
		os << pos;

		// Column 2: legacy-compatible total energies.
		bool wrote_total = false;
		for(size_t idx = 0; idx < n; idx++){
			const int len = access_lens[idx];
			if(len <= 0) continue;
			const int last_start = seq_len - len;
			if(pos > last_start) continue;
			const auto& vec = total_probs[idx];
			if(pos >= static_cast<int>(vec.size())) continue;
			const double energy = energy_of_prob(vec[pos], rt);
			os << (wrote_total ? "" : "\t") << len << "," << energy << ";";
			wrote_total = true;
		}

		// Column 3: by-loop energies (parser-safe extra column).
		bool wrote_byloop = false;
		for(size_t idx = 0; idx < n; idx++){
			const int len = access_lens[idx];
			if(len <= 0) continue;
			const int last_start = seq_len - len;
			if(pos > last_start) continue;
			if(pos >= static_cast<int>(exterior_probs[idx].size()) ||
			   pos >= static_cast<int>(hairpin_probs[idx].size()) ||
			   pos >= static_cast<int>(bulge_probs[idx].size()) ||
			   pos >= static_cast<int>(internal_probs[idx].size()) ||
			   pos >= static_cast<int>(multiloop_probs[idx].size())) {
				continue;
			}
		if(!wrote_byloop){
			os << "\tbyloop:";
			wrote_byloop = true;
		}
		os << len
		   << ",E=" << exterior_probs[idx][pos]
		   << ",H=" << hairpin_probs[idx][pos]
		   << ",B=" << bulge_probs[idx][pos]
		   << ",I=" << internal_probs[idx][pos]
		   << ",M=" << multiloop_probs[idx][pos]
		   << ";";
	}

		os << std::endl;
	}
	os << std::endl;
	os.precision(old_precision);
}

} // namespace io
} // namespace lcr

void LinCapR::output(std::ofstream& ofs, const std::string& seq_name) const {
	lcr::io::loop_profile_io(ofs, seq_name, prob_B, prob_E, prob_H, prob_I, prob_M, prob_S);
}
