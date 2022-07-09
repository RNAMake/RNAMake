//
// Created by Joe Yesselman on 6/28/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_STATE_BASEPAIR_HPP_
#define RNAMAKE_SRC_STRUCTURE_STATE_BASEPAIR_HPP_

#include <math/rotation.hpp>

namespace structure::state {
class Basepair {
public:
  Basepair(String &name, math::Vector3 &center, math::Vector3s &c1_prime_coords,
           math::Matrix3x3 &ref_frame)
      : _name(std::move(name)), _center(std::move(center)),
        _c1_prime_coords(std::move(c1_prime_coords)),
        _ref_frame(std::move(ref_frame)) {}

  Basepair(const Basepair &) = default;

  ~Basepair() = default;

public: // non const methods //////////////////////////////////////////////////
  void move(const math::Vector3 &p) {}


public:
  [[nodiscard]] inline const String &get_name() const { return _name; }

  [[nodiscard]] inline const math::Vector3 &get_center() const {
    return _center;
  }

  [[nodiscard]] inline const math::Vector3s &get_c1_prime_coords() const {
    return _c1_prime_coords;
  }

  [[nodiscard]] inline const math::Matrix3x3 &get_ref_frame() const {
    return _ref_frame;
  }

private:
  String _name;
  math::Vector3 _center;
  math::Vector3s _c1_prime_coords;
  math::Matrix3x3 _ref_frame;
};

typedef std::vector<Basepair> Basepairs;

} // namespace structure::state

#endif // RNAMAKE_SRC_STRUCTURE_STATE_BASEPAIR_HPP_
