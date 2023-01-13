//
// Created by Joseph Yesselman on 12/15/17.
//

#ifndef RNAMAKE_NEW_ALL_ATOM_BASEPAIR_H
#define RNAMAKE_NEW_ALL_ATOM_BASEPAIR_H

//#include <math/numerical.h>

#include <math/matrix_3x3.hpp>
#include <math/vector_3.hpp>
#include <structure/base/base.hpp>
#include <util/x3dna/x3dna.h>

namespace structure::all_atom {

class Basepair {
public:
  /// @brief - constructor
  inline Basepair(const util::Uuid &res1_uuid, util::Uuid const &res2_uuid,
                  const util::Uuid &uuid, structure::base::BasepairType bp_type,
                  util::x3dna::X3dnaBPType x3dna_type, String &name,
                  const math::Vector3 &center, math::Vector3s &c1_prime_coords,
                  const math::Matrix3x3 &ref_frame)
      : _res1_uuid(res1_uuid), _res2_uuid(res2_uuid), _uuid(uuid),
        _bp_type(bp_type), _x3dna_type(x3dna_type), _name(std::move(name)),
        _center(center), _c1_prime_coords(std::move(c1_prime_coords)),
        _ref_frame(ref_frame) {}

  Basepair(const Basepair &bp) = default;

  /// @brief - deconstructor
  ~Basepair() = default;

public:
  /// @brief - checks if two basepairs are equal
  inline bool is_equal(Basepair const &bp, bool check_uuid = true) const {
    if (check_uuid) {
      if (_uuid != bp._uuid) {
        return false;
      };
      if (_res1_uuid != bp._res1_uuid) {
        return false;
      }
      if (_res2_uuid != bp._res2_uuid) {
        return false;
      }
    }

    if (!math::are_points_equal(_center, bp._center)) {
      return false;
    }
    if (!math::are_matrices_equal(_ref_frame, bp._ref_frame)) {
      return false;
    }
    if (!math::are_points_equal(_c1_prime_coords[0], bp._c1_prime_coords[0])) {
      return false;
    }
    if (!math::are_points_equal(_c1_prime_coords[1], bp._c1_prime_coords[1])) {
      return false;
    }
    return true;
  }

public: // non const methods //////////////////////////////////////////////////
  /// @brief - moves the basepair by the specified vector
  void move(const math::Vector3 &p) {
    _center = _center + p;
    _c1_prime_coords[0] = _c1_prime_coords[0] + p;
    _c1_prime_coords[1] = _c1_prime_coords[1] + p;
  }
  /// @brief - rotates the basepair by the specified rotation matrix
  void rotate(const math::Matrix3x3 &rot) {
    _center = rot.dot(_center);
    _c1_prime_coords[0] = rot.dot(_c1_prime_coords[0]);
    _c1_prime_coords[1] = rot.dot(_c1_prime_coords[1]);
    _ref_frame = _ref_frame * rot.get_transposed();
    _ref_frame.unitarize();
  }
  /// @brief - swaps the residues
  inline void swap_residue_positions() {
    std::swap(_res1_uuid, _res2_uuid);
    std::swap(_c1_prime_coords[0], _c1_prime_coords[1]);
  }
  /// @brief - flips the reference frame
  inline void invert_reference_frame() {
    _ref_frame = _ref_frame.get_flip_orientation();
  }
  /// @brief - creates new uuids
  inline void new_uuids(const util::Uuid &r1_uuid, const util::Uuid &r2_uuid) {
    _res1_uuid = r1_uuid;
    _res2_uuid = r2_uuid;
    _uuid = util::generate_uuid();
  }

public:
  String get_str() const {
    String str = "";
    str += _center.get_str() +
           " " +
           _ref_frame.get_str() +
           " " +
           _c1_prime_coords[0].get_str() +
           " " +
           _c1_prime_coords[1].get_str() +
           " " +
           structure::base::BasepairTypeStrings[(int)_bp_type] +
           " " +
           std::to_string((int)_x3dna_type) +
           " " +
           _name;
    return str;
  }

public: // getters
  /// @brief - returns the uuid of the opposite base in the pair
  [[nodiscard]] const util::Uuid &get_partner(const util::Uuid &uuid) const {
    if (uuid == _res1_uuid) {
      return _res2_uuid;
    } else {
      return _res1_uuid;
    }
  }

  /// @brief - returns basepair type
  [[nodiscard]] inline structure::base::BasepairType get_bp_type() const {
    return _bp_type;
  }

  /// @brief - gets basepair uuid
  [[nodiscard]] inline const util::Uuid &get_uuid() const { return _uuid; }

  /// @brief - gets basepair name
  [[nodiscard]] inline const String &get_name() const { return _name; }

  /// @brief - gets the UUID of the 1st base in the pair
  [[nodiscard]] inline const util::Uuid &get_res1_uuid() const {
    return _res1_uuid;
  }

  /// @brief - gets the UUID of the 2nd base in the pair
  [[nodiscard]] inline util::Uuid const &get_res2_uuid() const {
    return _res2_uuid;
  }

  /// @brief - gets the reference frame
  [[nodiscard]] inline const math::Matrix3x3 &get_ref_frame() const {
    return _ref_frame;
  }

  /// @brief - gets the center of the basepair
  [[nodiscard]] inline const math::Vector3 &get_center() const {
    return _center;
  }

  /// @brief - gets the set of coordinates of the Carbon 1 primes
  [[nodiscard]] inline const math::Vector3s &get_c1_prime_coords() const {
    return _c1_prime_coords;
  }

  /// @brief - gets the coordinates of the carbon 1 prime of the first residue
  /// in the basepair
  [[nodiscard]] inline const math::Vector3 &get_res1_c1_prime_coord() const {
    return _c1_prime_coords[0];
  }

  /// @brief - gets the coordinates of the carbon 1 prime of the second residue
  /// in the basepair
  [[nodiscard]] inline const math::Vector3 &get_res2_c1_prime_coord() const {
    return _c1_prime_coords[1];
  }

private:
  util::Uuid _uuid;
  util::Uuid _res1_uuid;
  util::Uuid _res2_uuid;
  structure::base::BasepairType _bp_type;
  // TODO did we want to take this out? it's left out for now
  util::x3dna::X3dnaBPType _x3dna_type;
  String _name;
  math::Vector3 _center;
  math::Vector3s _c1_prime_coords;
  math::Matrix3x3 _ref_frame;
};

typedef std::shared_ptr<Basepair> BasepairOP;
typedef std::vector<Basepair> Basepairs;
typedef std::vector<BasepairOP> BasepairOPs;

/*
inline String generate_bp_name(Residue const &res1, Residue const &res2) {
  return primitives::generate_bp_name<Residue>(res1, res2);
}

primitives::BasepairType generate_bp_type(Residue const &, Residue const &,
                                          util::X3dnaBPType);  */

} // namespace structure::all_atom

#endif // RNAMAKE_NEW_ALL_ATOM_BASEPAIR_H