
//
//  residue.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_residue__
#define __RNAMake__ss_residue__

#include <cassert>
#include <cstdio>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <base/string.hpp>
#include <base/types.hpp>
#include <math/rotation.hpp>
#include <structure/base/base.hpp>

namespace structure::secondary_structure {

class Residue {
public:
  inline Residue(char name, char structure_code, int num, String &chain_id,
                 char const &i_code, util::Uuid const &uuid,
                 structure::base::ResidueType rtype)
      : _name(name), _structure_code(structure_code), _num(num),
        _chain_id(chain_id), _i_code(i_code), _uuid(uuid), _rtype(rtype) {
    _res_code = _assign_res_code(_name);
  }

  virtual ~Residue() = default;

public:
  /// @brief - checks if residues are equal
  inline bool operator==(const Residue &r) const { return is_equal(r); }

  /// @brief - checks if residues are not equal
  inline bool operator!=(const Residue &r) const { return !is_equal(r); }

public:
  /// @brief - checks if residues are equal (with slightly different arguments)
  [[nodiscard]] inline bool is_equal(const Residue &r,
                                     bool check_uuid = true) const {
    if (check_uuid && _uuid != r._uuid) {
      return false;
    }
    if (_name != r._name) {
      return false;
    }
    if (_structure_code != r._structure_code) {
      return false;
    }
    if (_num != r._num) {
      return false;
    }
    if (_chain_id != r._chain_id) {
      return false;
    }
    if (_i_code != r._i_code) {
      return false;
    }
    return true;
  }

public:
  /// @brief - stringifies the residue
  [[nodiscard]] inline String get_str() const {
    std::stringstream ss;
    // ss << _name << "," << dot_bracket_ << "," << _num << "," << _chain_id <<
    //     ","
    //    << _i_code;
    return ss.str();
  }

public: // getters

  /// @brief - gets the structure code
  [[nodiscard]] inline char get_structure_code() const {
    return _structure_code;
  }

  /// @brief - gets the residue code
  [[nodiscard]] inline int get_res_code() const { return _res_code; }

public: // setters

  /// @brief - sets the name of the residue code
  inline void set_name(char name) {
    _name = name;
    _res_code = _assign_res_code(_name);
  }

private:
  /// @brief - assigns number codes to eacah base in the residue
  int _assign_res_code(char name) {
    if (_rtype != structure::base::ResidueType::RNA) {
      return 99;
    }
    if (name == 'A') {
      return 0;
    } else if (name == 'C') {
      return 1;
    } else if (name == 'G') {
      return 2;
    } else if (name == 'U') {
      return 3;
    } else if (name == 'T') {
      return 3;
    } else if (name == 'N') {
      return -1;
    } else {
      throw structure::base::StructureException(
          "in sstruct::Residue encountered a unknown name: " +
          std::string(1, name));
    }
  }

private:
  // A = 0, C = 1, G = 2, U = 3, T = 3
  char _name;
  int _num;
  String _chain_id;
  char _i_code;
  util::Uuid _uuid;
  structure::base::ResidueType _rtype;
  char _structure_code;
  int _res_code;
};

/*struct res_less_than_key {
    inline
    bool
    operator() (ResidueOP const & r1, ResidueOP const & r2) {
        return (r1->num() < r2->num());
    }
};*/

typedef std::vector<Residue> Residues;

} // namespace structure::secondary_structure

#endif /* defined(__RNAMake__ss_residue__) */