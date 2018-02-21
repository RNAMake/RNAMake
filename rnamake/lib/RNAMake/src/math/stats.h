//
// Created by Joseph Yesselman on 2/20/18.
//

#ifndef TEST_STATS_H
#define TEST_STATS_H

#include <base/types.h>

double
sum(
        std::vector<double> const &);

double
sqsum(
        std::vector<double> const &);

double
stdev(
        std::vector<double> const &);

double
mean(
        std::vector<double> const &);

double
pearson_coeff(
        std::vector<double> const &,
        std::vector<double> const &);

double
avg_unsigned_diff(
        std::vector<double> const &,
        std::vector<double> const &);

#endif //TEST_STATS_H
