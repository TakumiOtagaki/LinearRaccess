#include <linearraccess/api.hpp>

#include <cmath>

int main() {
  const auto result = lcr::api::linear_raccess("GCGAAACGC", {1, 3});
  return result.accessibility.size() == 2
          && std::isfinite(result.accessibility[1][0])
      ? 0
      : 1;
}
