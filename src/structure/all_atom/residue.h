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

/// @brief - calculates the coordinates of the center of the residue
math::Vector3 center_of_atoms(const Atoms &);

class Residue {
public: // construction ///////////////////////////////////////////////////////
  /// @brief - constructors
  Residue(char name, int num, String &chain_id, char i_code, Atoms &atoms,
          const util::Uuid &uuid, structure::base::ResidueType rtype)
      : _name(name), _num(num), _chain_id(std::move(chain_id)), _i_code(i_code),
        _atoms(std::move(atoms)), _uuid(uuid), _rtype(rtype) {
    _build_beads();
  }

  Residue(const Residue &) = default;

  /// @brief - destructor
  ~Residue() = default;

public: // iterator ////////////////////////////////////////////////////////////
  typedef Atoms::const_iterator const_iterator;

  /// @brief - retrieves the beginning atom (?)
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
  /// @brief - are residues equal? (==)
  inline bool operator==(const Residue &r) const { return is_equal(r); }

  /// @brief - are residues not equal? (!=)
  inline bool operator!=(const Residue &r) const { return !is_equal(r); }

  friend std::ostream &operator<<(std::ostream &stream, const Residue &r) {
    /*stream << r._num << r._chain_id;
    if (r._i_code != ' ') {
      stream << "(" << r._i_code << ")";
    } */
    return stream;
  }

public:
  /// @brief - checks if two residues' components are equal
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
  /// @brief - moves residue by specified position vector
  void move(const math::Vector3 &p) {
    for (auto &a : _atoms) {
      a.move(p);
    }
    for (auto &b : _beads) {
      b.move(p);
    }
  }

  /// @brief - rotates residue by specified rotation matrix
  void rotate(const math::Matrix3x3 &rot) {
    for (auto &a : _atoms) {
      a.rotate(rot);
    }
    for (auto &b : _beads) {
      b.rotate(rot);
    }
  }

  /// @brief -
  inline void remove_beads() { _beads = util::Beads(); }

  /// @brief -
  inline void build_beads() { _build_beads(); }

  /// @brief - generates new UUID
  inline void new_uuid() { _uuid = util::generate_uuid(); }

public: // getters ////////////////////////////////////////////////////////////
  /// @brief - gets the name of the residue
  [[nodiscard]] inline char get_name() const { return _name; }

  /// @brief - gets the num of the residue
  [[nodiscard]] inline int get_num() const { return _num; }

  /// @brief - gets the id of the chain the residue is a part of
  [[nodiscard]] inline const String &get_chain_id() const { return _chain_id; }

  /// @brief - gets the i-code
  // TODO what's an i-code?
  [[nodiscard]] inline char get_i_code() const { return _i_code; }

  /// @brief - gets the UUID of the residue
  [[nodiscard]] inline util::Uuid get_uuid() const { return _uuid; }

  /// @brief - gets the residue type
  [[nodiscard]] inline structure::base::ResidueType get_rtype() const {
    return _rtype;
  }

  /// @brief - gets the beads
  [[nodiscard]] inline const util::Beads &get_beads() const { return _beads; }

  /// @brief - gets the specified atom in the residue
  [[nodiscard]] inline const Atom &get_atom(const String &name) const {
    for (auto const &a : _atoms) {
      if (a.get_name() == name) {
        return a;
      }
    }

    throw structure::base::StructureException(
        "atom name: " + name + " does not exist in this residue");
  }

  /// @brief - gets the coordinates of the specified atom
  [[nodiscard]] inline const math::Vector3 &
  get_coords(const String &name) const {
    return get_atom(name).get_coords();
  }

  // get coords of residue
  /// @brief - gets the coordinates of the center of the residue
  [[nodiscard]] math::Vector3 get_center() const;

  /// @brief - gets the x-coordinate of the center of the residue
  [[nodiscard]] double get_center_x() const;

  /// @brief - gets the y-coordinate of the center of the residue
  [[nodiscard]] double get_center_y() const;

  /// @brief - gets the z-coordinate of the center of the residue
  [[nodiscard]] double get_center_z() const;

  /// @brief - gets the number of atoms in the residue
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