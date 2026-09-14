#include "linearraccess/api.hpp"

#include "beam_prune.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void check_close(double actual, double expected, double tolerance,
                 const std::string& label) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(label + ": got " + std::to_string(actual)
                             + ", expected " + std::to_string(expected));
  }
}

lcr::api::LinearRaccessConfig exact_config(lcr::api::EnergyModel model,
                                            int sequence_length) {
  lcr::api::LinearRaccessConfig config;
  config.beam = 0;
  config.c_multi = sequence_length;
  config.c_hairpin = sequence_length;
  config.normalize_profiles = false;
  config.fast_logsumexp = false;
  config.energy_model = model;
  config.engine = lcr::api::EnergyEngine::LinearCapR;
  return config;
}

void check_model(lcr::api::EnergyModel model,
                 const std::vector<double>& golden_w1,
                 const std::vector<double>& golden_w3) {
  const std::string sequence = "GGGGGGAAAAAACCCCCC";
  std::vector<int> lengths{1, 2, 3, 4, 5};
  const auto config = exact_config(model, static_cast<int>(sequence.size()));
  const auto result = lcr::api::linear_raccess(sequence, lengths, config);
  const auto profile = lcr::api::lincapr_profile(sequence, config);

  require(result.sequence == sequence, "sequence normalization changed RNA input");
  require(result.accessibility.size() == lengths.size(), "length count mismatch");
  for (std::size_t length_index = 0; length_index < lengths.size(); ++length_index) {
    const int length = lengths[length_index];
    const std::size_t expected_size = sequence.size() - length + 1;
    require(result.accessibility[length_index].size() == expected_size,
            "window count mismatch");
    for (std::size_t i = 0; i < expected_size; ++i) {
      const double total = result.accessibility[length_index][i];
      const double components = result.exterior[length_index][i]
          + result.hairpin[length_index][i]
          + result.bulge[length_index][i]
          + result.internal[length_index][i]
          + result.multiloop[length_index][i];
      require(std::isfinite(total), "non-finite accessibility");
      require(total >= -1e-12 && total <= 1.0 + 1e-12,
              "probability outside numerical tolerance");
      check_close(total, components, 1e-12, "context decomposition");
      if (length_index + 1 < lengths.size()
          && i < result.accessibility[length_index + 1].size()) {
        const double longer = result.accessibility[length_index + 1][i];
        const std::string location = " model="
            + std::to_string(static_cast<int>(model))
            + " w=" + std::to_string(length)
            + " start=" + std::to_string(i);
        require(longer <= total + 1e-10,
                "left nesting monotonicity failed:" + location
                    + " longer=" + std::to_string(longer)
                    + " shorter=" + std::to_string(total)
                    + " excess=" + std::to_string(longer - total));
        const double right = result.accessibility[length_index][i + 1];
        require(longer <= right + 1e-10,
                "right nesting monotonicity failed:" + location
                    + " longer=" + std::to_string(longer)
                    + " shorter=" + std::to_string(right)
                    + " excess=" + std::to_string(longer - right));
      }
    }
  }

  for (std::size_t i = 0; i < sequence.size(); ++i) {
    check_close(result.accessibility[0][i], 1.0 - profile.stem[i], 1e-12,
                "P(w=1) identity");
  }

  const std::vector<std::size_t> starts{0, 5, 12, 17};
  for (std::size_t i = 0; i < starts.size(); ++i) {
    check_close(result.accessibility[0][starts[i]], golden_w1[i], 1e-12,
                "frozen w=1 value");
  }
  const std::vector<std::size_t> starts_w3{0, 5, 10, 15};
  for (std::size_t i = 0; i < starts_w3.size(); ++i) {
    check_close(result.accessibility[2][starts_w3[i]], golden_w3[i], 1e-12,
                "frozen w=3 value");
  }
}

void check_input_contract() {
  const auto config = exact_config(lcr::api::EnergyModel::Turner2004, 6);
  const auto canonical = lcr::api::linear_raccess("CCUAGG", {1}, config);
  const auto alias = lcr::api::linear_raccess("cctagg", {1}, config);
  require(alias.sequence == "CCUAGG", "DNA/lowercase normalization failed");
  for (std::size_t i = 0; i < canonical.accessibility[0].size(); ++i) {
    check_close(alias.accessibility[0][i], canonical.accessibility[0][i],
                1e-12, "normalized-input equivalence");
  }

  bool invalid_base_rejected = false;
  bool duplicate_length_rejected = false;
  bool negative_beam_rejected = false;
  try {
    static_cast<void>(lcr::api::linear_raccess("ACGN", {1}, config));
  } catch (const std::invalid_argument&) {
    invalid_base_rejected = true;
  }
  try {
    static_cast<void>(lcr::api::linear_raccess("ACGU", {1, 1}, config));
  } catch (const std::invalid_argument&) {
    duplicate_length_rejected = true;
  }
  try {
    auto invalid = config;
    invalid.beam = -1;
    static_cast<void>(lcr::api::linear_raccess("ACGU", {1}, invalid));
  } catch (const std::invalid_argument&) {
    negative_beam_rejected = true;
  }
  require(invalid_base_rejected && duplicate_length_rejected
              && negative_beam_rejected,
          "invalid API input was accepted");

#ifndef LINEARRACCESS_WITH_RACCESS
  bool unavailable_backend_rejected = false;
  try {
    auto unavailable = config;
    unavailable.engine = lcr::api::EnergyEngine::Raccess;
    static_cast<void>(lcr::api::linear_raccess("ACGU", {1}, unavailable));
  } catch (const std::invalid_argument&) {
    unavailable_backend_rejected = true;
  }
  require(unavailable_backend_rejected, "unavailable Raccess backend was accepted");
#endif
}

void check_beam_ties() {
  Map<int, Float> states;
  for (int i = 0; i < 5; ++i) states[i] = 0.0;
  lcr::beam::prune_states(states, 2,
                          [](int, Float score) { return score; });
  require(states.size() == 2 && states.count(0) == 1 && states.count(1) == 1,
          "beam tie-breaking is not deterministic");
}

void check_hairpin_cap_consistency() {
  const std::string sequence =
      "GUGCCCCCUCCGAUGCCGAGAGAUCCAAAGUCAUCGUUCCAAUGGGGGCUUUGUAGUCUG";
  auto config = exact_config(lcr::api::EnergyModel::Turner2004,
                             static_cast<int>(sequence.size()));
  config.c_hairpin = 10;
  const std::vector<int> lengths{1, 3, 7, 10, 20};
  const auto result = lcr::api::linear_raccess(sequence, lengths, config);
  const auto profile = lcr::api::lincapr_profile(sequence, config);
  for (std::size_t i = 0; i < sequence.size(); ++i) {
    check_close(profile.bulge[i], result.bulge[0][i], 1e-12,
                "hairpin cap must not omit bulge profile mass");
    check_close(profile.internal[i], result.internal[0][i], 1e-12,
                "hairpin cap must not omit internal profile mass");
    const double sum = profile.bulge[i] + profile.internal[i]
        + profile.hairpin[i] + profile.exterior[i]
        + profile.multiloop[i] + profile.stem[i];
    check_close(sum, 1.0, 1e-10, "raw capped structural profile normalization");
  }
  for (std::size_t length_index = 0; length_index < lengths.size(); ++length_index) {
    for (std::size_t i = 0; i < result.accessibility[length_index].size(); ++i) {
      const double probability = result.accessibility[length_index][i];
      require(probability >= -1e-12 && probability <= 1.0 + 1e-12,
              "hairpin cap produced an invalid probability");
      if (length_index + 1 < lengths.size()
          && i < result.accessibility[length_index + 1].size()) {
        const double longer = result.accessibility[length_index + 1][i];
        require(longer <= probability + 1e-10,
                "hairpin cap broke nested-window monotonicity");
      }
    }
  }
}

}  // namespace

int main() {
  try {
    check_model(
        lcr::api::EnergyModel::Turner2004,
        {0.025681506131667933, 0.065953600070891835,
         0.010443553608759194, 0.080893647475351299},
        {1.6951130550548978e-06, 0.065953600070891835,
         0.010443553608759194, 7.0185099888237464e-06});
    check_model(
        lcr::api::EnergyModel::Turner1999,
        {0.020028243222185535, 0.11234286570435763,
         0.060053828910026037, 0.072303475539784903},
        {1.1723816019069093e-06, 0.11234286570435763,
         0.060053828910026037, 4.194054355636042e-06});
    check_input_contract();
    check_beam_ties();
    check_hairpin_cap_consistency();
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  std::cout << "core_test\tpass\n";
  return EXIT_SUCCESS;
}
