//
// Created by Joseph Yesselman on 11/12/17.
//

#ifndef RNAMAKE_NEW_BEAD_H
#define RNAMAKE_NEW_BEAD_H

#include <base/string.hpp>
#include <math/matrix_3x3.hpp>
#include <math/vector_3.hpp>

namespace util {

/**
 * BeadType is an ENUM type. This is to specify which center of atoms each bead
 * represents.
 *
 * Phosphate (0):  P, OP1, OP2\n
 * Sugar (1):  O5',C5',C4',O4',C3',O3',C1',C2',O2'\n
 * Base  (2):  All remaining atoms
 */
enum class BeadType {
  PHOS,
  SUGAR,
  BASE,
  CALPHA,
  MCENTER // center of small molecule
};

/**
 * Exception for beads
 */
class BeadException : public std::runtime_error {
public:
  /**
   * Standard constructor for BeadException
   * @param   message   Error message for bead
   */
  explicit BeadException(String const &message) : std::runtime_error(message) {}
};

/**
 * Bead class stores information related to keeping track of steric clashes
 * between residues during building. They are never used outside the Residue
 * class
 *
 */
class Bead {
public:
  /**
   * Standard constructor for Bead object.
   * @param   btype   type of bead (PHOS, SUGAR or BASE)
   * @param   center  the average 3D position of all atoms represented by bead
   */
  inline Bead(const math::Vector3 &center, const BeadType bead_type)
      : _center(center), _bead_type(bead_type) {}

  inline Bead(const String &s) {
    auto spl = base::string::split(s, ",");
    _center = math::vector_from_str(spl[0]);
    _bead_type = BeadType(std::stoi(spl[1]));
  }

  inline Bead(const Bead &b) = default;

  ~Bead() = default;

public:
  [[nodiscard]] inline double distance(const Bead &b) const {
    return b._center.distance(_center);
  }

  inline void move(const math::Vector3 &p) { _center = _center + p; }

  inline void rotate(const math::Matrix3x3 &rot) { _center = rot.dot(_center); }

public: // getters
  [[nodiscard]] inline const math::Vector3 &get_center() const {
    return _center;
  }

  [[nodiscard]] inline BeadType get_type() const { return _bead_type; }

  [[nodiscard]] String get_str() const {
    return _center.get_str() + "," + std::to_string((int)_bead_type);
  }

  String get_type_name() {
    if (_bead_type == BeadType::PHOS) {
      return "PHOSPHATE";
    } else if (_bead_type == BeadType::SUGAR) {
      return "SUGAR";
    } else if (_bead_type == BeadType::BASE) {
      return "BASE";
    } else {
      throw BeadException("unknown bead type");
    }
  }

private:
  math::Vector3 _center;

  BeadType _bead_type;
};

typedef std::vector<Bead> Beads;

} // namespace util

#endif // RNAMAKE_NEW_BEAD_H