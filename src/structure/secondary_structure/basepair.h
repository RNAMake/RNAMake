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

  String get_str(Residue &res1, Residue &res2) const {
    // The way the old code wrote this is hidden in
    // secondary_structure::Motif::to_str
    String s;
    // Get residue positions
    // return "res1.position res2.position"
  }

  String name(Residue &res1, Residue &res2) const {
    std::stringstream ss;
    ss << res1.get_chain_id() << res1.get_num() << res1.get_i_code();
    String str1 = ss.str();
    std::stringstream ss2;
    ss2 << res2.get_chain_id() << res2.get_num() << res2.get_i_code();
    String str2 = ss2.str();
    if (str1 < str2) {
      return str1 + "-" + str2;
    } else {
      return str2 + "-" + str1;
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