#include "linearraccess/api.hpp"

#include "beam_inside_outside.hpp"
#ifdef LINEARRACCESS_WITH_RACCESS
#include "energy_raccess.hpp"
#endif
#include "linearraccess/dp_table_api.hpp"
#include "linearraccess/seq_utils.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lcr {
namespace api {

namespace {

energy::Model to_energy_model(EnergyModel model) {
	switch(model){
	case EnergyModel::Turner1999:
		return energy::Model::Turner1999;
	case EnergyModel::Turner2004:
		return energy::Model::Turner2004;
	}
	throw std::invalid_argument("unknown energy model");
}

LinCapR::EnergyEngine to_energy_engine(EnergyEngine engine) {
	switch(engine){
	case EnergyEngine::LinearCapR:
		return LinCapR::EnergyEngine::LinearCapR;
	case EnergyEngine::Raccess:
#ifdef LINEARRACCESS_WITH_RACCESS
		return LinCapR::EnergyEngine::Raccess;
#else
		throw std::invalid_argument(
		    "Raccess backend is unavailable in this build; rebuild with "
		    "LINEARRACCESS_WITH_RACCESS=ON and an authorized Raccess source tree");
#endif
	}
	throw std::invalid_argument("unknown energy engine");
}

double length_factor_for(EnergyEngine engine) {
#ifdef LINEARRACCESS_WITH_RACCESS
	return engine == EnergyEngine::Raccess ? raccess_length_factor() : 0.0;
#else
	static_cast<void>(engine);
	return 0.0;
#endif
}

void validate_config(const LinearRaccessConfig& config) {
	if(config.beam < 0 || config.c_multi < 0 || config.c_hairpin < 0){
		throw std::invalid_argument("beam, c_multi, and c_hairpin must be non-negative");
	}
	if(!(config.rt > 0.0) || !std::isfinite(config.rt)){
		throw std::invalid_argument("rt must be a finite positive output-conversion constant");
	}
	static_cast<void>(to_energy_model(config.energy_model));
	static_cast<void>(to_energy_engine(config.engine));
	if(config.raccess_use_lincapr_external
	   && config.engine != EnergyEngine::Raccess){
		throw std::invalid_argument(
		    "raccess_use_lincapr_external requires the Raccess engine");
	}
}

std::vector<int> validate_access_lens(const std::vector<int>& raw, int seq_n) {
	std::vector<int> out;
	out.reserve(raw.size());
	for(const int len : raw){
		if(len <= 0 || len > seq_n){
			throw std::invalid_argument(
			    "access length must be in [1, sequence length]: "
			    + std::to_string(len));
		}
		if(std::find(out.begin(), out.end(), len) != out.end()){
			throw std::invalid_argument(
			    "access lengths must be unique: " + std::to_string(len));
		}
		out.push_back(len);
	}
	return out;
}

} // namespace

LinearRaccessResult linear_raccess(const std::string& sequence,
                                   const std::vector<int>& access_lens,
                                   const LinearRaccessConfig& config) {
	validate_config(config);
	const std::string normalized_sequence = lcr::seq::normalize_sequence(sequence);
	LinearRaccessResult result;
	result.sequence = normalized_sequence;
	result.rt = config.rt;
	result.length_factor = length_factor_for(config.engine);
	const int seq_n = static_cast<int>(normalized_sequence.size());

	const auto lens = validate_access_lens(access_lens, seq_n);
	result.access_lens = lens;
	if(lens.empty()){
		return result;
	}

	LinCapR lcr(config.beam,
	            to_energy_model(config.energy_model),
	            to_energy_engine(config.engine),
	            config.raccess_use_lincapr_external,
	            config.normalize_profiles,
	            0.01,
	            config.c_multi,
	            config.c_hairpin);
	if(config.fast_logsumexp){
		lcr::dp::set_logsumexp_fast_mode();
	}else{
		lcr::dp::set_logsumexp_legacy_mode();
	}
	lcr.run(normalized_sequence);

	result.accessibility.reserve(lens.size());
	result.exterior.reserve(lens.size());
	result.hairpin.reserve(lens.size());
	result.bulge.reserve(lens.size());
	result.internal.reserve(lens.size());
	result.multiloop.reserve(lens.size());
	for(const int len : lens){
		const auto by = lcr.calc_accessibility_by_loop(len);
		result.accessibility.push_back(by.total);
		result.exterior.push_back(by.exterior);
		result.hairpin.push_back(by.hairpin);
		result.bulge.push_back(by.bulge);
		result.internal.push_back(by.internal);
		result.multiloop.push_back(by.multiloop);
	}

	return result;
}

LinCapRProfileResult lincapr_profile(const std::string& sequence,
                                     const LinearRaccessConfig& config) {
	validate_config(config);
	const std::string normalized_sequence = lcr::seq::normalize_sequence(sequence);
	LinCapRProfileResult result;
	result.sequence = normalized_sequence;

	LinCapR lcr(config.beam,
	            to_energy_model(config.energy_model),
	            to_energy_engine(config.engine),
	            config.raccess_use_lincapr_external,
	            config.normalize_profiles,
	            0.01,
	            config.c_multi,
	            config.c_hairpin);
	if(config.fast_logsumexp){
		lcr::dp::set_logsumexp_fast_mode();
	}else{
		lcr::dp::set_logsumexp_legacy_mode();
	}
	lcr.run(normalized_sequence);

	result.kT = lcr.get_kT();
	result.logZ = lcr.get_logZ();
	result.length_factor = length_factor_for(config.engine);
	result.bulge = lcr.get_prob_bulge();
	result.exterior = lcr.get_prob_exterior();
	result.hairpin = lcr.get_prob_hairpin();
	result.internal = lcr.get_prob_internal();
	result.multiloop = lcr.get_prob_multiloop();
	result.stem = lcr.get_prob_stem();
	return result;
}

} // namespace api
} // namespace lcr
