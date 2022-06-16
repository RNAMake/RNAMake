//
// Created by Joseph Yesselman on 2/20/18.
//

#include <cmath>

#include <math/stats.hpp>

namespace math {

double sum(std::vector<double> const &a) {
  double s = 0;
  for (int i = 0; i < a.size(); i++) {
    s += a[i];
  }
  return s;
}

double sqsum(std::vector<double> const &a) {
  double s = 0;
  for (int i = 0; i < a.size(); i++) {
    s += pow(a[i], 2);
  }
  return s;
}

double stdev(std::vector<double> const &nums) {
  double N = nums.size();
  return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

double mean(std::vector<double> const &a) { return sum(a) / a.size(); }

double pearson_coeff(std::vector<double> const &x,
                     std::vector<double> const &y) {
  double sum = 0;
  auto mean_x = mean(x);
  auto mean_y = mean(y);
  for (int i = 0; i < y.size(); i++) {
    sum += (x[i] - mean_x) * (y[i] - mean_y);
  }

  return sum / (x.size() * stdev(x) * stdev(y));

  // return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}

double avg_unsigned_diff(std::vector<double> const &x,
                         std::vector<double> const &y) {
  double diff = 0;
  for (int i = 0; i < x.size(); i++) {
    diff += std::abs(x[i] - y[i]);
  }
  return diff / x.size();
}

} // namespace math