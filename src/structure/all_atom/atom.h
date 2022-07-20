#ifndef __RNAMake__atom__
#define __RNAMake__atom__

#include <cstdio>

// RNAMake Headers
#include <base/exception.hpp>
#include <base/types.hpp>
#include <math/rotation.hpp>
#include <math/numerical.hpp>
#include <util/uuid.h>
#include <utility>

namespace structure::all_atom {
class Atom {
public:
  /// @brief - constructor
  inline Atom(String &name, math::Vector3 const &coords)
      : _name(std::move(name)), _coords(coords) {




  }

  inline explicit Atom(String const &s) {
    auto spl = ::base::string::split(s, " ");
    if (spl.size() != 4) {
      throw ::base::InputException("tried to initialize atom with string: " +
                                   s);
    }
    _name = spl[0];
    _coords =
        math::Vector3(std::stod(spl[1]), std::stod(spl[2]), std::stod(spl[3]));
  }
  /**
   * Copy constructor
   * @param   a   atom object to from
   */
  inline Atom(Atom const &a) = default;

public:
  inline bool operator==(Atom const &a) const {
    if (_name != a._name) {
      return false;
    }
    if (!math::are_points_equal(_coords, a._coords)) {
       return false;
    }
    return true;
  }

  inline bool operator!=(Atom const &a) const { return !(*this == a); }

public: // non const methods //////////////////////////////////////////////////
  /// @brief - moves atoms by distance "p"; the distance is added to coords
  inline void move(const math::Vector3 &p) { _coords = _coords + p; }

  /// @brief - rotates a point by a degree defined by a rotation matrix
  inline void rotate(const math::Matrix3x3 &rot) {
    _coords = rot.dot(_coords);
  }

public: // trival getters /////////////////////////////////////////////////////
  /// @brief - gets the name of the atom
  [[nodiscard]] inline const String &get_name() const { return _name; }

  /// @brief - gets the coords of the atom
  [[nodiscard]] inline const math::Vector3 &get_coords() const {
    return _coords;
  }

public: // coord getters //////////////////////////////////////////////////////
  /// @brief - gets x coordinate of the atom
  [[nodiscard]] inline double get_x() const { return _coords.get_x(); }
  /// @brief - gets y coordinate of the atom
  [[nodiscard]] inline double get_y() const { return _coords.get_y(); }
  /// @brief - gets z coordinate of the atom
  [[nodiscard]] inline double get_z() const { return _coords.get_z(); }

private:
  String _name;
  math::Vector3 _coords;
};

// TODO Need to remove these pointers
//  (do we really?)
typedef std::vector<Atom> Atoms;

} // namespace structure::all_atom

#endif /* defined(__RNAMake__atom__) */