#include "beam_inside_outside.hpp"
#ifdef LINEARRACCESS_WITH_RACCESS
#include "energy_raccess.hpp"
#endif
#include "FileReader.hpp"
#include "io.hpp"
#include "linearraccess/dp_table_api.hpp"
#include "linearraccess/version.hpp"
#include "build_revision.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

constexpr double kDefaultOutputRtKcalMol = 0.61633008;

bool parse_integer(const std::string& raw, int& value){
	try{
		size_t used = 0;
		const long parsed = std::stol(raw, &used);
		if(used != raw.size()
		   || parsed < std::numeric_limits<int>::min()
		   || parsed > std::numeric_limits<int>::max()) return false;
		value = static_cast<int>(parsed);
		return true;
	}catch(const std::exception&){
		return false;
	}
}

bool parse_real(const std::string& raw, double& value){
	try{
		size_t used = 0;
		value = std::stod(raw, &used);
		return used == raw.size();
	}catch(const std::exception&){
		return false;
	}
}

bool parse_access_lens(const std::string& raw, std::vector<int>& values){
	values.clear();
	std::string tmp = raw;
	for(char& c : tmp){
		if(c == ';' || c == ':' || c == ',') c = ' ';
	}
	std::istringstream iss(tmp);
	std::string token;
	while(iss >> token){
		int val = 0;
		if(!parse_integer(token, val) || val <= 0) return false;
		if(std::find(values.begin(), values.end(), val) != values.end()) return false;
		values.push_back(val);
	}
	return !values.empty();
}

bool is_metadata_token(const std::string& value){
	return !value.empty() && std::all_of(value.begin(), value.end(), [](char c) {
		const unsigned char byte = static_cast<unsigned char>(c);
		return std::isalnum(byte) || c == '.' || c == '_' || c == '-' || c == '+';
	});
}

const char* model_name(energy::Model model){
	switch(model){
	case energy::Model::Turner1999:
		return "turner1999";
	case energy::Model::Turner2004:
		return "turner2004";
	}
	throw std::invalid_argument("unknown energy model");
}

const char* engine_name(LinCapR::EnergyEngine engine){
	switch(engine){
	case LinCapR::EnergyEngine::LinearCapR:
		return "lincapr";
	case LinCapR::EnergyEngine::Raccess:
		return "raccess";
	}
	throw std::invalid_argument("unknown energy engine");
}

double engine_length_factor(LinCapR::EnergyEngine engine){
#ifdef LINEARRACCESS_WITH_RACCESS
	if(engine == LinCapR::EnergyEngine::Raccess){
		return lcr::raccess_length_factor();
	}
#else
	static_cast<void>(engine);
#endif
	return 0.0;
}

double folding_rt_kcal_mol(energy::Model model,
							LinCapR::EnergyEngine engine){
#ifdef LINEARRACCESS_WITH_RACCESS
	if(engine == LinCapR::EnergyEngine::Raccess){
		return lcr::raccess_rt_kcal_mol();
	}
#else
	static_cast<void>(engine);
#endif
	// LinearCapR stores energies and kT in 10 cal/mol units.
	return energy::get_params(model).kT / 100.0;
}

void write_metadata(std::ostream& os,
					const std::string& source_revision,
					const std::vector<int>& access_lens,
					int beam,
					int c_multi,
					int c_hairpin,
					double output_rt,
					bool normalize_profiles,
					bool fast_logsumexp,
					bool output_probabilities,
					bool byloop,
					bool hybrid,
					energy::Model model,
					LinCapR::EnergyEngine engine){
	const std::streamsize old_precision = os.precision();
	os << std::setprecision(std::numeric_limits<double>::max_digits10);
	os << "# linearraccess_metadata=1\n";
	os << "# software_version=" << lcr::kSoftwareVersion << "\n";
	os << "# source_revision=" << source_revision << "\n";
	os << "# engine=" << engine_name(engine) << "\n";
	os << "# requested_energy_model=" << model_name(model) << "\n";
	os << "# effective_energy_model="
	   << (engine == LinCapR::EnergyEngine::Raccess ? "raccess-turner1999"
	                                               : model_name(model))
	   << "\n";
	os << "# beam=" << beam << "\n";
	os << "# c_multi=" << c_multi << "\n";
	os << "# c_hairpin=" << c_hairpin << "\n";
	os << "# logsumexp=" << (fast_logsumexp ? "fast" : "exact") << "\n";
	os << "# normalize_profiles=" << (normalize_profiles ? "true" : "false") << "\n";
	os << "# length_factor=" << engine_length_factor(engine) << "\n";
	os << "# folding_rt=" << folding_rt_kcal_mol(model, engine) << "\n";
	os << "# output_rt=" << output_rt << "\n";
	os << "# output_mode="
	   << (output_probabilities ? "probability" : (byloop ? "energy-by-loop" : "energy"))
	   << "\n";
	os << "# experimental_hybrid=" << (hybrid ? "true" : "false") << "\n";
	os << "# access_lens=";
	for(size_t i = 0; i < access_lens.size(); ++i){
		if(i > 0) os << ',';
		os << access_lens[i];
	}
	os << "\n";
	os.precision(old_precision);
}

void usage(){
	std::cout << "Usage: ./LinRacc -seqfile=<fa> -outfile=<txt> -access_len=1,5,10 [options]\n";
	std::cout << "Options:\n";
	std::cout << "  -beam=<int>        Beam size (default: 200)\n";
	std::cout << "  -energy=<model>    lincapr engine only: turner2004 (default) or turner1999\n";
	std::cout << "  -engine=<name>     lincapr (default)";
#ifdef LINEARRACCESS_WITH_RACCESS
	std::cout << " or raccess (optional; fixed Turner 1999)";
#endif
	std::cout << "\n";
	std::cout << "  -c_multi=<int>     Max multiloop unpaired length (default: 30)\n";
	std::cout << "  -c_hairpin=<int>   Max direct hairpin-loop length (default: 30)\n";
	std::cout << "  -no-normalize      Disable per-position profile normalization\n";
	std::cout << "  -no-fast-logsumexp Disable polynomial logsumexp approximation\n";
	std::cout << "  -rt=<float>        RT used only for probability-to-energy output (default: 0.61633008)\n";
	std::cout << "  -probabilities     Output raw probabilities instead of -RT log(P)\n";
	std::cout << "  -byloop            Output per-loop accessibility breakdown (extra column)\n";
	std::cout << "  -metadata          Prepend reproducibility metadata comments\n";
	std::cout << "  -source-revision=<id>  Revision token recorded with -metadata\n";
	std::cout << "  --build-info          Print embedded build version/revision as JSON\n";
	std::cout << "  -debug-exit        Call std::exit(0) after run (skip destructors)\n";
	std::cout << "  -debug-release-energy  Reset energy model after run\n";
	std::cout << "  -debug-leak-energy     Leak energy model after run\n";
	std::cout << "  -debug-skip-access     Skip accessibility calc and exit after run\n";
	std::cout << "  -debug-outer-dp-now i0,i1  Dump alpha_O/beta_O right after run()\n";
	std::cout << "  -debug-beam j0,j1      Log beam stats for j in [j0,j1]\n";
	std::cout << "  -debug-beam-interval n Log beam stats every n positions\n";
	std::cout << "  -debug-beam-stop       Stop after j1 when debug-beam is set\n";
	std::cout << "  -debug-access          Log accessibility calc progress\n";
#ifdef LINEARRACCESS_WITH_RACCESS
	std::cout << "  -raccess-use-lincapr-external EXPERIMENTAL hybrid; do not use for reference results\n";
#endif
}

} // namespace

int main(int argc, char** argv){
	if(argc == 2 && std::string(argv[1]) == "--build-info"){
		std::cout << "{\"software_version\":\"" << lcr::kSoftwareVersion
		          << "\",\"source_revision\":\"" << lcr::kBuildSourceRevision << "\"}\n";
		return 0;
	}
	if(argc < 2){
		usage();
		return 1;
	}

	std::string seqfile;
	std::string outfile;
	std::vector<int> access_lens;
	int beam = 200;
	int c_multi = 30;
	int c_hairpin = MAXLOOP;
	double rt = kDefaultOutputRtKcalMol;
	bool normalize_profiles = true;
	bool disable_fast_logsumexp = false;
	bool output_probabilities = false;
	bool byloop = false;
	bool output_metadata = false;
	std::string source_revision = lcr::kBuildSourceRevision;
	bool debug_exit = false;
	bool debug_release_energy = false;
	bool debug_leak_energy = false;
	bool debug_skip_access = false;
	bool debug_outer_dp_now = false;
	int debug_outer_dp_now_i0 = -1;
	int debug_outer_dp_now_i1 = -1;
	bool debug_beam = false;
	int debug_beam_j0 = -1;
	int debug_beam_j1 = -1;
	int debug_beam_interval = 1;
	bool debug_beam_stop = false;
	bool debug_access = false;
	bool raccess_use_lincapr_external = false;
	energy::Model energy_model = energy::Model::Turner2004;
	LinCapR::EnergyEngine engine = LinCapR::EnergyEngine::LinearCapR;

	for(int i = 1; i < argc; i++){
		std::string arg = argv[i];
		if(arg.rfind("-seqfile=", 0) == 0){
			seqfile = arg.substr(strlen("-seqfile="));
		}else if(arg.rfind("-outfile=", 0) == 0){
			outfile = arg.substr(strlen("-outfile="));
		}else if(arg.rfind("-access_len=", 0) == 0){
			if(!parse_access_lens(arg.substr(strlen("-access_len=")), access_lens)){
				std::cout << "Error: -access_len must be a list of positive integers.\n";
				return 1;
			}
		}else if(arg.rfind("-beam=", 0) == 0){
			if(!parse_integer(arg.substr(strlen("-beam=")), beam)){
				std::cout << "Error: -beam must be an integer.\n";
				return 1;
			}
		}else if(arg.rfind("-c_multi=", 0) == 0){
			if(!parse_integer(arg.substr(strlen("-c_multi=")), c_multi)){
				std::cout << "Error: -c_multi must be an integer.\n";
				return 1;
			}
		}else if(arg.rfind("-c_hairpin=", 0) == 0){
			if(!parse_integer(arg.substr(strlen("-c_hairpin=")), c_hairpin)){
				std::cout << "Error: -c_hairpin must be an integer.\n";
				return 1;
			}
		}else if(arg.rfind("-rt=", 0) == 0){
			if(!parse_real(arg.substr(strlen("-rt=")), rt)){
				std::cout << "Error: -rt must be a number.\n";
				return 1;
			}
		}else if(arg == "-no-normalize"){
			normalize_profiles = false;
		}else if(arg == "-no-fast-logsumexp"){
			disable_fast_logsumexp = true;
		}else if(arg == "-probabilities" || arg == "--probabilities"){
			output_probabilities = true;
		}else if(arg == "-byloop" || arg == "--byloop"){
			byloop = true;
		}else if(arg == "-metadata" || arg == "--metadata"){
			output_metadata = true;
		}else if(arg.rfind("-source-revision=", 0) == 0
		         || arg.rfind("--source-revision=", 0) == 0){
			const size_t equals = arg.find('=');
			source_revision = arg.substr(equals + 1);
			if(!is_metadata_token(source_revision)){
				std::cout << "Error: source revision must be a non-empty identifier token.\n";
				return 1;
			}
			output_metadata = true;
		}else if(arg == "-debug-exit"){
			debug_exit = true;
		}else if(arg == "-debug-release-energy"){
			debug_release_energy = true;
		}else if(arg == "-debug-leak-energy"){
			debug_leak_energy = true;
		}else if(arg == "-debug-skip-access"){
			debug_skip_access = true;
		}else if(arg == "-debug-outer-dp-now"){
			if(i + 1 >= argc){
				std::cout << "Error: -debug-outer-dp-now requires i0,i1" << std::endl;
				return 1;
			}
			const std::string val = argv[++i];
			const size_t comma = val.find(',');
			if(comma == std::string::npos){
				std::cout << "Error: -debug-outer-dp-now expects i0,i1" << std::endl;
				return 1;
			}
			debug_outer_dp_now_i0 = std::stoi(val.substr(0, comma));
			debug_outer_dp_now_i1 = std::stoi(val.substr(comma + 1));
			debug_outer_dp_now = true;
		}else if(arg == "-debug-beam"){
			if(i + 1 >= argc){
				std::cout << "Error: -debug-beam requires j0,j1" << std::endl;
				return 1;
			}
			const std::string val = argv[++i];
			const size_t comma = val.find(',');
			if(comma == std::string::npos){
				std::cout << "Error: -debug-beam expects j0,j1" << std::endl;
				return 1;
			}
			debug_beam_j0 = std::stoi(val.substr(0, comma));
			debug_beam_j1 = std::stoi(val.substr(comma + 1));
			debug_beam = true;
		}else if(arg == "-debug-beam-interval"){
			if(i + 1 >= argc){
				std::cout << "Error: -debug-beam-interval requires n" << std::endl;
				return 1;
			}
			debug_beam_interval = std::max(1, std::atoi(argv[++i]));
		}else if(arg == "-debug-beam-stop"){
			debug_beam_stop = true;
		}else if(arg == "-debug-access"){
			debug_access = true;
		}else if(arg == "-raccess-use-lincapr-external" || arg == "--raccess-use-lincapr-external"){
#ifdef LINEARRACCESS_WITH_RACCESS
			raccess_use_lincapr_external = true;
#else
			std::cout << "Error: Raccess backend is unavailable in this build.\n";
			return 1;
#endif
		}else if(arg.rfind("-energy=", 0) == 0){
			const std::string val = arg.substr(strlen("-energy="));
			if(val == "turner1999"){
				energy_model = energy::Model::Turner1999;
			}else if(val == "turner2004"){
				energy_model = energy::Model::Turner2004;
			}else{
				std::cout << "Error: invalid energy model: " << val << std::endl;
				return 1;
			}
		}else if(arg.rfind("-engine=", 0) == 0){
			const std::string val = arg.substr(strlen("-engine="));
			if(val == "lincapr"){
				engine = LinCapR::EnergyEngine::LinearCapR;
			}else if(val == "raccess"){
#ifdef LINEARRACCESS_WITH_RACCESS
				engine = LinCapR::EnergyEngine::Raccess;
#else
				std::cout << "Error: Raccess backend is unavailable in this build.\n";
				return 1;
#endif
			}else{
				std::cout << "Error: invalid energy engine: " << val << std::endl;
				return 1;
			}
		}else if(arg == "-h" || arg == "--help"){
			usage();
			return 0;
		}else{
			std::cout << "Unknown option: " << arg << std::endl;
			usage();
			return 1;
		}
	}

	if(seqfile.empty() || outfile.empty() || access_lens.empty()){
		std::cout << "Error: -seqfile, -outfile, -access_len are required.\n";
		usage();
		return 1;
	}
	if(beam < 0 || c_multi < 0 || c_hairpin < 0){
		std::cout << "Error: -beam, -c_multi, and -c_hairpin must be non-negative.\n";
		return 1;
	}
	if(!(rt > 0.0) || !std::isfinite(rt)){
		std::cout << "Error: -rt must be finite and positive.\n";
		return 1;
	}
	if(output_probabilities && byloop){
		std::cout << "Error: -probabilities and -byloop cannot be combined.\n";
		return 1;
	}
	if(raccess_use_lincapr_external && engine != LinCapR::EnergyEngine::Raccess){
		std::cout << "Error: -raccess-use-lincapr-external requires -engine=raccess.\n";
		return 1;
	}

	FileReader reader;
	std::vector<std::string> seqs;
	std::vector<std::string> names;
	if(!reader.read(seqfile, seqs, names)){
		return 1;
	}
	for(size_t idx = 0; idx < seqs.size(); ++idx){
		for(const int len : access_lens){
			if(len > static_cast<int>(seqs[idx].size())){
				std::cout << "Error: access length " << len
				          << " exceeds sequence length for " << names[idx] << ".\n";
				return 1;
			}
		}
	}

	std::ofstream ofs(outfile);
	if(!ofs){
		std::cout << "Error: cannot open output file: " << outfile << std::endl;
		return 1;
	}
	if(output_metadata){
		write_metadata(ofs, source_revision, access_lens, beam, c_multi, c_hairpin,
		               rt, normalize_profiles, !disable_fast_logsumexp,
		               output_probabilities, byloop, raccess_use_lincapr_external,
		               energy_model, engine);
	}

	for(size_t idx = 0; idx < seqs.size(); idx++){
		LinCapR lcr(beam, energy_model, engine, raccess_use_lincapr_external,
		            normalize_profiles, 0.01, c_multi, c_hairpin);
		if(disable_fast_logsumexp){
			lcr::dp::set_logsumexp_legacy_mode();
		}
		if(debug_beam){
			lcr.set_debug_beam(debug_beam_j0, debug_beam_j1, debug_beam_interval, debug_beam_stop);
		}
		lcr.run(seqs[idx]);
		if(debug_outer_dp_now){
			lcr.debug_outer_dp_dump(debug_outer_dp_now_i0, debug_outer_dp_now_i1);
		}
		if(debug_release_energy){
			lcr.release_energy(false);
		}
		if(debug_leak_energy){
			lcr.release_energy(true);
		}
		if(debug_exit){
			std::exit(0);
		}
		if(debug_skip_access){
			return 0;
		}

		std::vector<std::vector<double>> total_probs;
		total_probs.reserve(access_lens.size());

		std::vector<std::vector<double>> exterior_probs;
		std::vector<std::vector<double>> hairpin_probs;
		std::vector<std::vector<double>> bulge_probs;
		std::vector<std::vector<double>> internal_probs;
		std::vector<std::vector<double>> multiloop_probs;
		if(byloop){
			exterior_probs.reserve(access_lens.size());
			hairpin_probs.reserve(access_lens.size());
			bulge_probs.reserve(access_lens.size());
			internal_probs.reserve(access_lens.size());
			multiloop_probs.reserve(access_lens.size());
		}

		for(const int len : access_lens){
			if(debug_access){
				std::cerr << "debug_access: len=" << len << " start" << std::endl;
			}
			if(byloop){
				const auto by = lcr.calc_accessibility_by_loop(len);
				if(debug_access){
					std::cerr << "debug_access: len=" << len
					          << " total_size=" << by.total.size()
					          << " exterior_size=" << by.exterior.size()
					          << std::endl;
				}
				total_probs.push_back(by.total);
				exterior_probs.push_back(by.exterior);
				hairpin_probs.push_back(by.hairpin);
				bulge_probs.push_back(by.bulge);
				internal_probs.push_back(by.internal);
				multiloop_probs.push_back(by.multiloop);
			}else{
				auto probs = lcr.calc_accessibility(len);
				if(debug_access){
					std::cerr << "debug_access: len=" << len
					          << " total_size=" << probs.size()
					          << std::endl;
				}
				total_probs.push_back(std::move(probs));
			}
			if(debug_access){
				std::cerr << "debug_access: len=" << len << " done" << std::endl;
			}
		}

		if(byloop){
			lcr::io::accessibility_raccess_byloop_io(
				ofs,
				names[idx],
				static_cast<int>(seqs[idx].size()),
				access_lens,
				total_probs,
				exterior_probs,
				hairpin_probs,
				bulge_probs,
				internal_probs,
				multiloop_probs,
				rt);
		}else if(output_probabilities){
			lcr::io::accessibility_probability_io(
				ofs,
				names[idx],
				static_cast<int>(seqs[idx].size()),
				access_lens,
				total_probs);
		}else{
			lcr::io::accessibility_raccess_io(
				ofs,
				names[idx],
				static_cast<int>(seqs[idx].size()),
				access_lens,
				total_probs,
				rt);
		}
	}

	return 0;
}
