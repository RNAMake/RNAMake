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
#include <base/types.hpp>
#include <structure/all_atom/atom.h>
#include <structure/base.hpp>
#include <util/bead.h>

namespace structure::all_atom {

math::Vector3 center(Atoms const &);

class Residue {
public:
  /**
   * Standard constructor for Residue object.
   * @param   name    name of residue i.e. "G", "C", etc
   * @param   num     residue num
   * @param   chain_id    the id of the chain i.e. "A", "B", only one character
   * @param   i_code  residue insertion code, usually nothing ("")
   * @param   rtype   residue type, dictates residue topology and information
   *
   */
  Residue(char name, int num, String &chain_id, char i_code, Atoms &atoms,
          const util::Uuid &uuid)
      : _name(name), _num(num), _chain_id(std::move(chain_id)), _i_code(i_code),
        _atoms(std::move(atoms)), _uuid(uuid) {
    //_build_beads();
  }

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
  inline bool operator==(Residue const &r) const { return is_equal(r); }

  inline bool operator!=(Residue const &r) const { return !is_equal(r); }

  friend std::ostream &operator<<(std::ostream &stream, Residue const &r) {
    /*stream << r._num << r._chain_id;
    if (r._i_code != ' ') {
      stream << "(" << r._i_code << ")";
    } */
    return stream;
  }
  // TODO I think this can be externalized
  inline int connected_to(Residue const &res, float cutoff = 5.0) const {
    String o3 = "O3'", p = "P";

    // 5' to 3'
    Atom o3_atom = get_atom(o3), p_atom = res.get_atom(p);
    if (o3_atom.get_coords().distance(p_atom.get_coords()) < cutoff) {
      return 1;
    }

    // 3' to 5'
    o3_atom = res.get_atom(o3);
    p_atom = get_atom(p);
    if (o3_atom.get_coords().distance(p_atom.get_coords()) < cutoff) {
      return -1;
    }

    return 0;
  }

public:
  inline bool is_equal(Residue const &r, bool check_uuid = true) const {
    /*if (check_uuid && _uuid != r._uuid) {
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
    }  */
    return true;
  }

public: // getters ////////////////////////////////////////////////////////////
  [[nodiscard]] char get_name() const { return _name; }

  [[nodiscard]] int get_num() const { return _num; }

  [[nodiscard]] const String &get_chain_id() const { return _chain_id; }

  [[nodiscard]] char get_i_code() const { return _i_code; }

  [[nodiscard]] const util::Uuid get_uuid() const { return _uuid; }

  [[nodiscard]] inline const Atom &get_atom(String const &name) const {
    for (auto const &a : _atoms) {
      if (a.get_name() == name) {
        return a;
      }
    }
    throw StructureException("atom name: " + name +
                             " does not exist in this residue");
  }

  [[nodiscard]] inline const Atom &get_atom(Index index) const {
    return _atoms[index];
  }

  [[nodiscard]] inline const math::Vector3 &
  get_coords(String const &name) const {
    return get_atom(name).get_coords();
  }

  [[nodiscard]] inline const util::Bead &
  get_bead(util::BeadType bead_type) const {
    for (auto const &b : _beads) {
      if (b.get_type() == bead_type) {
        return b;
      }
    }
    // throw ResidueException("bead type: " + std::to_string((int)bead_type) +
    //                        " does not exist in this residue");
    return _beads[0];
  }

  [[nodiscard]] inline math::Vector3 get_center() const {
    return center(_atoms);
  }

  [[nodiscard]] inline size_t get_num_atoms() const { return _atoms.size(); }

  // inline String const &get_res_name() const { return _res_type.get_name(); }

  // inline SetType get_res_set_type() const { return _res_type.get_set_type();
  // }

  String get_str() const;



  /**
   * wrapper for to_pdb_str(int &, int, String const &) when one does not care
   *about renumbering atoms and residue
   **/
  inline String get_pdb_str(int acount) const {
    // auto num = _num;
    // auto chain_id = _chain_id;
    // return get_pdb_str(acount, num, chain_id);
    return String("");
  }

  /**
   * returns pdb formatted string of residue's coordinate information
   * @param   acount  current atom index, default: 1
   * @param   rnum    starting residue number
   * @param   chain_id    the chain id of the chain, i.e. "A", "B" etc
   * @code
   *  #include "instances/structure_instances.hpp"
   *  auto r = instances::residue();
   *  //get PDB formatted String (can also write to file with
   * r->to_pdb("test.pdb") ) std::cout << r->to_pdb_str() << std::endl; //OUTPUT
   * --> ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00 ATOM
   * 2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00 ATOM      3 C4'  G
   * A 103     -24.521 -48.156  86.068  1.00  0.00 ATOM      4 O4'  G   A 103
   * -24.861 -49.568  86.118  1.00  0.00 ATOM      5 C3'  G   A 103     -23.009
   * -48.119  86.281  1.00  0.00 ATOM      6 O3'  G   A 103     -22.548
   * -46.872  86.808  1.00  0.00
   *  .
   *  .
   *  .
   * @endcode
   */
  String get_pdb_str(int &, int, char) const;

  String get_bead_pdb_str(int &, int, char) const;

  /**
   * writes a PDB string formmated verision of this Residue object to file
   * @param  filename of output PDB file
   */
  void write_pdb(String const);

public: // non const
  /*inline void move(math::Vector3 const &p) {
    for (auto &a : _atoms) {
      a.move(p);
    }
    for (auto &b : _beads) {
      b.move(p);
    }
  }

  inline void transform(math::Matrix3x3 const &r, math::Vector3 const &t,
                        math::Vector3 &dummy) {
    for (auto &a : _atoms) {
      a.transform(r, t, dummy);
    }
    for (auto &b : _beads) {
      b.transform(r, t, dummy);
    }
  }

  inline void transform(math::Matrix3x3 const &r, math::Vector3 const &t) {
    auto dummy = math::Vector3();
    transform(r, t, dummy);
  }               */

  inline void remove_beads() { _beads = util::Beads(); }

  inline void build_beads() { _build_beads(); }

  // inline void new_uuid() { _uuid = util::Uuid(); }


private:
  void _build_beads();

  void _build_beads_RNA();

  // TODO Need to remove this function later

private:
  char _name;
  int _num;
  String _chain_id;
  char _i_code;
  util::Uuid _uuid;
  /// @brief vector of the atom objects that belong to this residue
  Atoms _atoms;
  /// @brief vector of bead objects for sterics
  util::Beads _beads;
};

/**
 * Shared pointer typedef for Residue. Only use shared pointers!
 */
typedef std::vector<Residue> Residues;
typedef std::shared_ptr<Residue> ResidueOP;
typedef std::shared_ptr<Residue const> ResidueCOP;

/**
 * Typedef of a vector of shared pointer vectors, only used this.
 */
inline bool residue_steric_clash_RNA(Residue const &r1, Residue const &r2) {
  for (auto it1 = r1.bead_begin(); it1 != r1.bead_end(); it1++) {
    if (it1->get_type() == util::BeadType::PHOS) {
      continue;
    }
    for (auto it2 = r2.bead_begin(); it2 != r2.bead_end(); it2++) {
      if (it2->get_type() == util::BeadType::PHOS) {
        continue;
      }
      if (it1->distance(*it2) < 2.5) {
        return true;
      }
    }
  }

  /*std::for_each(r1.bead_begin(), r1.bead_end(), [](util::Bead const & b1) {
      if(b1.get_type() != util::BeadType::PHOS) {
          std::for_each(r2.bead_begin(), r2.bead_end(), [](util::Bead const &
  b2) { if (b2.get_type() != util::BeadType::PHOS) { if (b1.distance(b2) < 2.5)
  { return true; }
              }
          });
      }
  });*/

  return false;
}

inline bool residue_steric_clash(Residue const &r1, Residue const &r2) {
  for (auto it1 = r1.bead_begin(); it1 != r1.bead_end(); it1++) {
    for (auto it2 = r2.bead_begin(); it2 != r2.bead_end(); it2++) {
      if (it1->distance(*it2) < 2.5) {
        return true;
      }
    }
  }
  return false;
}

//  ResidueOP
//  _setup_residue(
//          String const & key,
//          Atoms const & atoms,
//          ResidueTypeCOP res_type)  {
//
//
//      auto spl = base::split_str_by_delimiter(key, "|");
//      auto atom_ptrs = std::vector<Atom const *>(res_type->get_num_atoms());
//
//      for(auto const & a : atoms) {
//          // not a valid atom for this residue
//          if(! res_type->is_valid_atom_name(a.get_name())) {
////                  LOGW <<  a.get_name() + " does not belong to residue " +
/// res_type->get_name() + ": IGNORING!";
//              continue;
//          }
//          auto index = res_type->get_atom_index(a.get_name());
//          atom_ptrs[index] = &a;
//      }
//
//      if(res_type->get_set_type() == SetType::RNA) {
//          auto i = 1;
//          auto missing_phosphate = false;
//          for (auto const & a : atom_ptrs) {
//              i++;
//              if (i < 5 && a == nullptr) { missing_phosphate = 1; }
//          }
//
//          if (missing_phosphate) {
//              if (!_replace_missing_phosphate_backbone(atom_ptrs, res_type)) {
////                      LOGW << "tried to fill in missing phosphate backbone
/// for residue " + spl[0] + " " + spl[1];
//              }
//          }
//      }
} // namespace structure::all_atom

#endif /* defined(__RNAMake__structure_residue__) */