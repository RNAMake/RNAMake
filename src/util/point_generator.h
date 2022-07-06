//
// Created by Joe Yesselman on 6/29/22.
//

#ifndef RNAMAKE_POINT_GENERATOR_H
#define RNAMAKE_POINT_GENERATOR_H

#include <math/vector_3.hpp>
#include <util/random_number_generator.h>

namespace util {
class PointGenerator {
public:
  PointGenerator()
      : _test_p(math::Vector3()), _rng(util::RandomNumberGenerator()) {}

  ~PointGenerator() {}

public:
  inline math::Vector3 const &rand_point(int scale) {
    _test_p.set_x(_rng.rand() * _rng.randrange(scale));
    _test_p.set_y(_rng.rand() * _rng.randrange(scale));
    _test_p.set_z(_rng.rand() * _rng.randrange(scale));

    if (_rng.randrange(1000) < 500) {
      _test_p.set_x(-_test_p.get_x());
    }
    if (_rng.randrange(1000) < 500) {
      _test_p.set_y(-_test_p.get_y());
    }
    if (_rng.randrange(1000) < 500) {
      _test_p.set_z(-_test_p.get_z());
    }

    return _test_p;
  }

private:
  math::Vector3 _test_p;
  util::RandomNumberGenerator _rng;
};
}
#endif // RNAMAKE_POINT_GENERATOR_H
