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
#include <structure/base/base.hpp>
#include <util/bead.h>

namespace structure::all_atom {

math::Vector3 center_of_atoms(const Atoms &atoms) {
  auto center = math::Vector3();
  for (auto const &a : atoms) {
    center += a.get_coords();
  }
  return center / float(atoms.size());
}

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
  inline bool operator==(Residue const &r) const { return is_equal(r); }

  inline bool operator!=(Residue const &r) const { return !is_equal(r); }

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
      /*if (_atoms[i] != r._atoms[i]) {
        return false;
      } */
    }
    return true;
  }

public: // getters ////////////////////////////////////////////////////////////
  [[nodiscard]] char get_name() const { return _name; }

  [[nodiscard]] int get_num() const { return _num; }

  [[nodiscard]] const String &get_chain_id() const { return _chain_id; }

  [[nodiscard]] char get_i_code() const { return _i_code; }

  [[nodiscard]] util::Uuid get_uuid() const { return _uuid; }

  [[nodiscard]] structure::base::ResidueType get_rtype() const {
    return _rtype;
  }

  [[nodiscard]] inline const util::Beads &get_beads() const { return _beads; }

  [[nodiscard]] inline const Atom &get_atom(String const &name) const {
    for (auto const &a : _atoms) {
      if (a.get_name() == name) {
        return a;
      }
    }
    throw structure::base::StructureException("atom name: " + name +
                             " does not exist in this residue");
  }

  [[nodiscard]] inline const math::Vector3 &
  get_coords(const String &name) const {
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
    return center_of_atoms(_atoms);
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
  void _build_beads() {
    if (_rtype == structure::base::ResidueType::RNA) {
      _build_beads_RNA();
    } else if (_rtype == structure::base::ResidueType::PROTEIN) {
      _beads.push_back(util::Bead(get_coords("CA"), util::BeadType::CALPHA));
    } else {
      _beads.push_back(util::Bead(get_center(), util::BeadType::MCENTER));
    }
  }

  void _build_beads_RNA() {
    std::vector<Atom const *> phos_atoms, sugar_atoms, base_atoms;
    int i = -1;
    for (auto const &a : _atoms) {
      i++;
      if (a.get_name() == "P" || a.get_name() == "OP1" ||
          a.get_name() == "OP2") {
        phos_atoms.push_back(&a);
      } else if (a.get_name().find('\'') != -1) {
        sugar_atoms.push_back(&a);
      } else {
        base_atoms.push_back(&a);
      }
    }
    auto get_center = [&](std::vector<Atom const *> const &atom_ptrs) {
      auto center = math::Vector3();
      for (auto a : atom_ptrs) {
        center = center + a->get_coords();
      }
      center = center / float(atom_ptrs.size());
      return center;
    };
    if (!phos_atoms.empty()) {
      _beads.push_back(
          util::Bead(get_center(phos_atoms), util::BeadType::PHOS));
    }
    if (!sugar_atoms.empty()) {
      _beads.push_back(
          util::Bead(get_center(sugar_atoms), util::BeadType::SUGAR));
    }
    if (!base_atoms.empty()) {
      _beads.push_back(
          util::Bead(get_center(base_atoms), util::BeadType::BASE));
    }
  }

public:
  void move(const math::Vector3 &p) {
    for (auto &a : _atoms) {
      a.move(p);
    }
  }

  void transform(const math::RotandTrans &rt) {
    for (auto &a : _atoms) {
      a.transform(rt);
    }
  }

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
Residue get_residue_from_str(const String &s) {
  Strings spl = ::base::string::split(s, ",");
  int a = 0;
  char name = spl[1][0];
  int num = std::stoi(spl[2]);
  String chain_id = spl[3];
  char i_code = ' ';
  util::Uuid uuid = util::generate_uuid();
  Atoms atoms;
  int i = 5;
  while (i < spl.size()) {
    if (spl[i].size() <= 1) {
      i++;
      continue;
    }
    try {
      atoms.push_back(Atom(spl[i]));
    } catch (...) {
      std::cout << "FAILED!" << std::endl;
      std::cout << i << " " << spl[i] << std::endl;
      for (int j = 5; j < spl.size(); j++) {
        std::cout << spl[j].length() << " ";
      }
      std::cout << std::endl;
      exit(0);
    }
    i++;
  }
  return {name, num, chain_id, i_code, atoms, uuid,
          structure::base::ResidueType::RNA};
}

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

} // namespace structure::all_atom

#endif /* defined(__RNAMake__structure_residue__) */