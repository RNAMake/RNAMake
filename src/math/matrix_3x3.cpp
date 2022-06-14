//
// Created by Joe Yesselman on 6/7/22.
//

#include <string>


#include "matrix_3x3.hpp"
#include <base/exception.hpp>
#include <base/string.hpp>

namespace math {

Matrix3x3 matrix_from_str(const String &s) {
  Strings doubles = base::string::split(s, " ");
  Reals m = {};
  if (doubles.empty()) {
    String msg = "No input detected!";
    base::log_and_throw<base::InputException>(msg);
  }

  if (doubles.size() != 9) {
    String msg = "Requires nine input arguments. From string got " + std::to_string(doubles.size());
    base::log_and_throw<base::InputException>(msg);
  }

  for (auto &i : doubles) {m.push_back(std::stod(i));
  }
  // check its 9 numbers
  return Matrix3x3{m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]};
}

}