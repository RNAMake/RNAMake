//
//  residue.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__structure_residue__
#define __RNAMake__structure_residue__

#include <cstdio>

// RNAMake Headers
#include <base/global_constants.hpp>
#include <base/types.hpp>
#include <structure/all_atom/atom.h>
#include <structure/base/base.hpp>
#include <util/bead.h>

namespace structure::all_atom {

math::Vector3 center_of_atoms(const Atoms &);

class Residue {
public: // construction ///////////////////////////////////////////////////////
  Residue(char name, int num, String &chain_id, char i_code, Atoms &atoms,
          const util::Uuid &uuid, structure::base::ResidueType rtype)
      : _name(name), _num(num), _chain_id(std::move(chain_id)), _i_code(i_code),
        _atoms(std::move(atoms)), _uuid(uuid), _rtype(rtype) {
    _build_beads();
  }

  Residue(const Residue &) = default;

  ~Residue() = default;

public: // iterator ////////////////////////////////////////////////////////////
  typedef Atoms::const_iterator const_iterator;

  [[nodiscard]] const_iterator begin() const noexcept { return _atoms.begin(); }

  [[nodiscard]] const_iterator end() const noexcept { return _atoms.end(); }

  typedef util::Beads::const_iterator bead_const_iterator;

  [[nodiscard]] bead_const_iterator bead_begin() const noexcept {
    return _beads.begin();
  }

  [[nodiscard]] bead_const_iterator bead_end() const noexcept {
    return _beads.end();
  }

public:
  inline bool operator==(const Residue &r) const { return is_equal(r); }

  inline bool operator!=(const Residue &r) const { return !is_equal(r); }

  friend std::ostream &operator<<(std::ostream &stream, const Residue &r) {
    /*stream << r._num << r._chain_id;
    if (r._i_code != ' ') {
      stream << "(" << r._i_code << ")";
    } */
    return stream;
  }

public:
  [[nodiscard]] inline bool is_equal(const Residue &r,
                                     bool check_uuid = true) const {
    if (check_uuid && _uuid != r._uuid) {
      return false;
    }
    if (_name != r._name) {
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
    for (int i = 0; i < _atoms.size(); i++) {
      if (_atoms[i] != r._atoms[i]) {
        return false;
      }
    }
    return true;
  }

public: // non const methods /////////////////////////////////////////////////
  void move(const math::Vector3 &p) {
    for (auto &a : _atoms) {
      a.move(p);
    }
    for(auto &b : _beads) {
      b.move(p);
    }
  }

  void rotate(const math::Matrix3x3 &rot) {
    for (auto &a : _atoms) {
      a.rotate(rot);
    }
    for (auto &b : _beads) {
      b.rotate(rot);
    }
  }

  inline void remove_beads() { _beads = util::Beads(); }

  inline void build_beads() { _build_beads(); }

  inline void new_uuid() { _uuid = util::generate_uuid(); }

public: // getters ////////////////////////////////////////////////////////////
  [[nodiscard]] inline char get_name() const { return _name; }

  [[nodiscard]] inline int get_num() const { return _num; }

  [[nodiscard]] inline const String &get_chain_id() const { return _chain_id; }

  [[nodiscard]] inline char get_i_code() const { return _i_code; }

  [[nodiscard]] inline util::Uuid get_uuid() const { return _uuid; }

  [[nodiscard]] inline structure::base::ResidueType get_rtype() const {
    return _rtype;
  }

  [[nodiscard]] inline const util::Beads &get_beads() const { return _beads; }

  [[nodiscard]] inline const Atom &get_atom(const String &name) const {
    for (auto const &a : _atoms) {
      if (a.get_name() == name) {
        return a;
      }
    }
    throw structure::base::StructureException(
        "atom name: " + name + " does not exist in this residue");
  }

  [[nodiscard]] inline const math::Vector3 &
  get_coords(const String &name) const {
    return get_atom(name).get_coords();
  }

  [[nodiscard]] math::Vector3 get_center() const;

  [[nodiscard]] inline size_t get_num_atoms() const { return _atoms.size(); }

  // inline String const &get_res_name() const { return _res_type.get_name(); }

  // inline SetType get_res_set_type() const { return _res_type.get_set_type();
  // }

  // String get_str() const;

private:
  void _build_beads();

  void _build_beads_RNA();

private:
  char _name;
  int _num;
  String _chain_id;
  char _i_code;
  util::Uuid _uuid;
  /// @brief vector of the atom objects that belong to this residue
  Atoms _atoms;
  structure::base::ResidueType _rtype;
  /// @brief vector of bead objects for sterics
  util::Beads _beads = {};
};

/**
 * Shared pointer typedef for Residue. Only use shared pointers!
 */
typedef std::vector<Residue> Residues;

// TODO do not return RNA as type always!
Residue get_residue_from_str(const String &);

inline bool residue_steric_clash(const Residue &r1, const Residue &r2) {
  for (auto it1 = r1.bead_begin(); it1 != r1.bead_end(); it1++) {
    for (auto it2 = r2.bead_begin(); it2 != r2.bead_end(); it2++) {
      if (it1->distance(*it2) < STERIC_CLASH_RADIUS) {
        return true;
      }
    }
  }
  return false;
}

} // namespace structure::all_atom

#endif /* defined(__RNAMake__structure_residue__) */