#pragma once

#include <string>
#include <vector>

namespace lcr {
namespace api {

enum class EnergyModel {
	Turner2004,
	Turner1999,
};

enum class EnergyEngine {
	Raccess,
	LinearCapR,
};

struct LinearRaccessConfig {
	int beam = 100;
	int c_multi = 30;
	int c_hairpin = 30;
	double rt = 0.61633008;
	bool normalize_profiles = true;
	bool fast_logsumexp = true;
	bool raccess_use_lincapr_external = false;
	EnergyModel energy_model = EnergyModel::Turner2004;
	EnergyEngine engine = EnergyEngine::LinearCapR;
};

struct LinearRaccessResult {
	std::string sequence;
	std::vector<int> access_lens;
	std::vector<std::vector<double>> accessibility;
	std::vector<std::vector<double>> exterior;
	std::vector<std::vector<double>> hairpin;
	std::vector<std::vector<double>> bulge;
	std::vector<std::vector<double>> internal;
	std::vector<std::vector<double>> multiloop;
	double rt = 0.61633008;
	double length_factor = 0.0;
};

LinearRaccessResult linear_raccess(const std::string& sequence,
                                   const std::vector<int>& access_lens,
                                   const LinearRaccessConfig& config = {});

struct LinCapRProfileResult {
	std::string sequence;
	double kT = 0.0;
	double logZ = 0.0;
	double length_factor = 0.0;
	std::vector<double> bulge;
	std::vector<double> exterior;
	std::vector<double> hairpin;
	std::vector<double> internal;
	std::vector<double> multiloop;
	std::vector<double> stem;
};

LinCapRProfileResult lincapr_profile(const std::string& sequence,
                                     const LinearRaccessConfig& config = {});

} // namespace api
} // namespace lcr
