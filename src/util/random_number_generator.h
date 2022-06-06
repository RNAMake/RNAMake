
//
// Created by Joseph Yesselman on 12/17/17.
//

#ifndef RNAMAKE_NEW_RANDOM_NUMBER_GENERATOR_H
#define RNAMAKE_NEW_RANDOM_NUMBER_GENERATOR_H

#include <random>

#include <math/vector_3.hpp>

namespace util {

  class RandomNumberGenerator {
  public:
      RandomNumberGenerator() {
          srand(unsigned(time(NULL)));
          std::random_device rd;
          std::mt19937 mt(rd());
          std::uniform_real_distribution<double> dist(0, 1);

          mt_ = mt;
          dist_ = dist;
      }

      inline
      double
      rand() { return dist_(mt_); }

      inline
      int
      randrange(int i) { return (int) (i * rand()); }



  private:
      std::mt19937 mt_;
      std::uniform_real_distribution<double> dist_;

  };

//  math::Point
//  get_random_point(
//          RandomNumberGenerator & rng,
//          int bound) {
//      auto x = bound - rng.rand()*2*bound;
//      auto y = bound - rng.rand()*2*bound;
//      auto z = bound - rng.rand()*2*bound;
//      return math::Point(x, y, z);
//  }


}

#endif //RNAMAKE_NEW_RANDOM_NUMBER_GENERATOR_H