//
// Created by Joe Yesselman on 5/25/22.
//
// std headers
#include <string>

// RNAMake Headers
#include "vector_3.hpp"
#include <base/exception.hpp>
#include <base/string.hpp>

namespace math {
void vectors_from_str(const String &s, math::Vector3s &vecs /* return */) {
  Strings doubles = base::string::split(s, " ");
  Reals point = {0, 0, 0};
  int pos = 0;

  if (s.empty()) {
    String msg = "No vector input detected!";
    base::log_and_throw<base::InputException>(msg);
  }

  for (auto const &d : doubles) {
    point[pos] = std::stod(d);
    pos += 1;
    if (pos == 3) {
      vecs.push_back(math::Vector3{point});
      pos = 0;
    }
  }
  // there are leftovers!
  if (pos != 0) {
    String msg = "vector must contain 3 numbers. From string got " +
                 std::to_string(doubles.size());
    base::log_and_throw<base::InputException>(msg);
  }
}
} // namespace math
