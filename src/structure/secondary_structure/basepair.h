//
// Created by Joseph Yesselman on 10/26/17.
//

#ifndef RNAMAKE_NEW_BASEPAIR_H
#define RNAMAKE_NEW_BASEPAIR_H

#include <base/string.hpp>
#include <structure/base/base.hpp>
#include <structure/secondary_structure/residue.h>

namespace structure::secondary_structure {

class Basepair {
public:
  inline Basepair(const util::Uuid &res1_uuid, const util::Uuid &res2_uuid,
                  const util::Uuid &uuid,
                  const structure::base::BasepairType bp_type)
      : _res1_uuid(res1_uuid), _res2_uuid(res2_uuid), _uuid(uuid),
        _bp_type(bp_type) {}

  Basepair(const Basepair &) = default;

  ~Basepair() = default;

public: // operators //////////////////////////////////////////////////////////
  inline bool operator==(const Basepair &bp) const {
    if (_res1_uuid != bp._res1_uuid) {
      return false;
    }
    if (_res2_uuid != bp._res2_uuid) {
      return false;
    }
    if (_uuid != bp._uuid) {
      return false;
    }
    if (_bp_type != bp._bp_type) {
      return false;
    }
    return true;
  }

  inline bool operator!=(Basepair const &bp) const { return !(*this == bp); }

public:
  /// @brief - checks if two basepairs are equal
  inline bool is_equal(Basepair const &bp, bool check_uuid = true) const {
    if (check_uuid && _res1_uuid != bp._res1_uuid) {
      return false;
    }
    if (check_uuid && _res2_uuid != bp._res2_uuid) {
      return false;
    }
    if (check_uuid && _uuid != bp._uuid) {
      return false;
    }
    if (_bp_type != bp._bp_type) {
      return false;
    }
    return true;
  }

public: // non const methods //////////////////////////////////////////////////


public: // getters ////////////////////////////////////////////////////////////
  /// @brief - gets second base in base pair
  [[nodiscard]] util::Uuid const &get_partner(util::Uuid const &uuid) const {
    if (uuid == _res1_uuid) {
      return _res2_uuid;
    } else {
      return _res1_uuid;
    }
  }

public: // trivial getters ///////////////////////////////////////////////////
  /// @brief - gets the basepair type
  [[nodiscard]] inline structure::base::BasepairType const &
  get_bp_type() const {
    return _bp_type;
  }

  /// @brief - gets the UUID of the basepair
  [[nodiscard]] inline util::Uuid const &get_uuid() const { return _uuid; }

  /// @brief - gets the UUID of an individual base
  [[nodiscard]] inline util::Uuid const &get_res1_uuid() const {
    return _res1_uuid;
  }

  /// @brief - gets the UUID of an individual base
  [[nodiscard]] inline util::Uuid const &get_res2_uuid() const {
    return _res2_uuid;
  }

private:
  util::Uuid _uuid;
  util::Uuid _res1_uuid;
  util::Uuid _res2_uuid;
  structure::base::BasepairType _bp_type;
};

typedef std::vector<Basepair> Basepairs;

} // namespace structure::secondary_structure

#endif // RNAMAKE_NEW_BASEPAIR_H