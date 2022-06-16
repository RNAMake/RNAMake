//
//  residue.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__structure_residue__
#define __RNAMake__structure_residue__

#include <stdio.h>

// RNAMake Headers
#include <base/types.hpp>
#include <doctest.h>
#include <primitives/residue.h>
#include <structure/atom.h>
#include <structure/residue_type.h>
#include <structure/residue_type_set.h>
#include <util/bead.h>
#include <util/uuid.h>

namespace structure {

math::Vector3 center(Atoms const &);

/**
 * Store residue information from pdb file, stores all Atom objects that
 * belong to residue. Implementation is designed to be extremely lightweight.
 *
 * To get an example instance of Residue please include:\n
 * #include "instances/structure_instances.hpp"
 *
 * Example of Common Usage:
 *
 * @code
 *  //grabbing example instance
 *  #include "instances/structure_instances.hpp"
 *  auto r = instances::residue();
 *  //simply getting name of residue, etc: A,G,C or U
 *  std::cout << r->name() << std::endl; //OUTPUT: "G"
 *
 *  //getting a specific atom from residue
 *  auto a = r->gFet_atom("C1'");
 *  std::cout << a->coords() << std::endl; //OUTPUT: "-23.806 -50.289  86.732"
 *
 *  auto beads = r->get_beads();
 *  std::cout << beads[0]->btype() << " " << beads[0]->center() << std::endl;
 *  //OUTPUT: "1 24.027 -48.5001111111 86.368"
 *
 *  //get PDB formatted String (can also write to file with
 * r->to_pdb("test.pdb") ) std::cout << r->to_pdb_str() << std::endl; //OUTPUT
 * --> ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00 ATOM
 * 2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00 ATOM      3 C4'  G
 * A 103     -24.521 -48.156  86.068  1.00  0.00F ATOM      4 O4'  G   A 103
 * -24.861 -49.568  86.118  1.00  0.00 ATOM      5 C3'  G   A 103     -23.009
 * -48.119  86.281  1.00  0.00 ATOM      6 O3'  G   A 103     -22.548
 * -46.872  86.808  1.00  0.00
 *  .
 *  .
 *  .
 * @endcode
 */
class Residue : public primitives::Residue {
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
  Residue(char name, int num, char chain_id, char i_code,
          ResidueTypeCOP res_type, Atoms &atoms, util::Uuid const &uuid)
      : primitives::Residue(name, num, chain_id, i_code, uuid),
        _res_type(res_type), _atoms(std::move(atoms)) {
    _build_beads();
  }

  /**
   * Copy constructor
   * @param   r   Residue object copying from
   */
  Residue(Residue const &r)
      : primitives::Residue(r._name, r._num, r._chain_id, r._i_code, r._uuid),
        _res_type(r._res_type), _atoms(r._atoms), _beads(r._beads) {}

  /**
   * Construction from String, used in reading data from files
   * @param   s   string generated from to_str()
   * @see to_str()
   */
  Residue(String const &s, structure::ResidueTypeSet const &rts)
      : primitives::Residue() {

    auto spl = base::string::split(s, ",");
    _res_type = rts.get_residue_type(spl[0]);
    _name = spl[1][0];
    _num = std::stoi(spl[2]);
    _chain_id = spl[3][0];
    _i_code = spl[4][0];
    _uuid = util::Uuid();
    _atoms = Atoms();
    int i = 5;
    while (i < spl.size()) {
      if (spl[i].size() <= 1) {
        i++;
        continue;
      }
      _atoms.push_back(Atom(spl[i]));
      i++;
    }
    //          _replace_missing_phosphate_backbone(_atoms, rts);
    _build_beads();
  }

  /**
   * Empty deconstructor
   */
  ~Residue() {}

public: // iterator stuff
  typedef Atoms::const_iterator const_iterator;

  const_iterator begin() const noexcept { return _atoms.begin(); }

  const_iterator end() const noexcept { return _atoms.end(); }

  typedef util::Beads::const_iterator bead_const_iterator;

  bead_const_iterator bead_begin() const noexcept { return _beads.begin(); }

  bead_const_iterator bead_end() const noexcept { return _beads.end(); }

public:
  inline bool operator==(Residue const &r) const { return is_equal(r); }

  inline bool operator!=(Residue const &r) const { return !is_equal(r); }

  friend std::ostream &operator<<(std::ostream &stream, Residue const &r) {
    stream << r._num << r._chain_id;
    if (r._i_code != ' ') {
      stream << "(" << r._i_code << ")";
    }
    return stream;
  }

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

public: // getters
  /**
   * get atom object by its name
   * @param   name   name of atom
   *
   * examples:
   * @code
   *  auto a = r->get_atom("C1'");
   *  std::cout << a->coords() << std::endl;
   *  //OUTPUT -23.806 -50.289  86.732
   * @endcode
   */
  inline Atom const &get_atom(String const &name) const {
    int index = _res_type->get_atom_index(name);
    if (index == -1) {
      throw ResidueException("atom name: " + name +
                             " does not exist in this residue");
    }
    return _atoms[index];
  }

  inline Atom const &get_atom(Index index) const { return _atoms[index]; }

  inline math::Vector3 const &get_coords(String const &name) const {
    return get_atom(name).get_coords();
  }

  inline util::Bead const &get_bead(util::BeadType bead_type) const {
    for (auto const &b : _beads) {
      if (b.get_type() == bead_type) {
        return b;
      }
    }
    throw ResidueException("bead type: " + std::to_string((int)bead_type) +
                           " does not exist in this residue");
  }

  // TODO Temporary get rid of these setters

  template <typename T> inline void set_chain_id(T &id) { _chain_id = id; }

  inline void set_num(int &num) { _num = num; }

  inline void set_uuid(util::Uuid const &id) { _uuid = id; }

  // TODO need to remove this function
  inline AtomOPs get_atoms() {
    AtomOPs atom_ptrs;
    for (auto const &a : _atoms) {
      AtomOP atom_ptr = std::make_shared<Atom>(a);
      atom_ptrs.push_back(atom_ptr);
    }
    return atom_ptrs;
  }

  inline math::Vector3 get_center() const { return center(_atoms); }

  inline size_t get_num_atoms() const { return _atoms.size(); }

  inline String const &get_res_name() const { return _res_type->get_name(); }

  inline SetType get_res_set_type() const { return _res_type->get_set_type(); }

  /**
   * stringifes residue object
   * @code
   *  #include "instances/structure_instances.hpp"
   *  auto r = instances::residue();
   *  std::cout << r->to_str() << std::endl;
   *  //OUTPUT
   *  GUA,G,103,A,,N,N,N,O5' -26.469 -47.756 84.669,C5' -25.05
   * -47.579 84.775,C4' -24.521 -48.156 86.068,O4' -24.861 -49.568 86.118,C3'
   * -23.009 -48.119 86.281,O3' -22.548 -46.872 86.808,C1' -23.806
   * -50.289 86.732,C2' -22.812 -49.259 87.269,O2' -23.167 -48.903 88.592,N1
   * -19.538 -52.485 85.025,C2 -19.717 -51.643 86.097,N2 -18.624
   * -51.354 86.809,N3 -20.884 -51.124 86.445,C4 -21.881 -51.521 85.623,C5
   * -21.811 -52.356 84.527,C6 -20.546 -52.91 84.164,O6 -20.273
   * -53.677 83.228,N7 -23.063 -52.513 83.947,C8 -23.858 -51.786 84.686,N9
   * -23.21 -51.159 85.722,
   * @endcode
   */
  String get_str() const;

  //      json::JSON
  //      get_json() const;

  /**
   * wrapper for to_pdb_str(int &, int, String const &) when one does not care
   *about renumbering atoms and residue
   **/
  inline String get_pdb_str(int acount) const {
    auto num = _num;
    auto chain_id = _chain_id;
    return get_pdb_str(acount, num, chain_id);
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
  inline void move(math::Vector3 const &p) {
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
  }

  inline void remove_beads() { _beads = util::Beads(); }

  inline void build_beads() { _build_beads(); }

  inline void new_uuid() { _uuid = util::Uuid(); }

public: // getters
        /**
         * getter for the one letter residue type
         */
        /*inline
        String
        short_name() const { return res_type_->short_name(); }*/

private:
  void _build_beads();

  void _build_beads_RNA();

  // TODO Need to remove this function later

  std::map<String, std::shared_ptr<Residue>> _ref_residues;

  //      bool
  //      _replace_missing_phosphate_backbone(
  //              std::vector<Atom> &atoms,
  //              ResidueTypeCOP res_type) const {
  //
  //          auto ref_res = _ref_residues.at(res_type->get_name());
  //
  //          // if these atoms do not exist cannot build res ref frame
  ////          if(atoms[ res_type->get_atom_index("C1'")] == nullptr) { return
  ///false; } /          if (res_type->get_short_name() == 'A' ||
  ///res_type->get_short_name() == 'G') { /              if(atoms[
  ///res_type->get_atom_index("N9")] == nullptr) { return false; } /          }
  ////          else {
  ////              if(atoms[ res_type->get_atom_index("N1")] == nullptr) {
  ///return false; } /          }
  //            bool found_C1 = false;
  //          for (Atom const &atom : atoms) {
  //              if(atom.get_name() == "C1'") {
  //                  found_C1 = true;
  //              }
  //          }
  //
  //          auto ref_frame_1 = _get_res_ref_frame_from_atoms(atoms, res_type);
  //          auto ref_frame_2 = _get_res_ref_frame(ref_res);
  //          auto rot = dot(ref_frame_1.transpose(), ref_frame_2);
  //          auto r_t = rot.transpose();
  //          auto t = -ref_res->get_center();
  //          auto c4p_atom = atoms[res_type->get_atom_index("C4'")];
  //          if (c4p_atom == nullptr) { return false; }
  //
  //          ref_res->transform(r_t, t);
  //          ref_res->move(c4p_atom.get_coords() - ref_res->get_coords("C4'"));
  //
  //          for (int i = 0; i < 5; i++) {
  //              atoms[i] = new Atom(ref_res->get_atom(i).get_name(),
  //                                  ref_res->get_atom(i).get_coords());
  //          }
  //          return true;
  //      }
  //
  //      math::Matrix
  //      _get_res_ref_frame(std::shared_ptr<Residue const> r) const {
  //          auto vec1 = math::Point();
  //          auto vec2 = math::Point();
  //          if (r->get_name() == 'A' || r->get_name() == 'G') {
  //              vec1 = (r->get_coords("N9") -
  //              r->get_coords("C1'")).normalize(); vec2 = (r->get_coords("N9")
  //              - r->get_bead(util::BeadType::BASE).get_center()).normalize();
  //          } else {
  //              vec1 = (r->get_coords("N1") -
  //              r->get_coords("C1'")).normalize(); vec2 = (r->get_coords("N1")
  //              - r->get_bead(util::BeadType::BASE).get_center()).normalize();
  //          }
  //          auto cross = vec1.cross(vec2);
  //          auto m = math::Matrix(
  //                  vec1.get_x(), vec1.get_y(), vec1.get_z(),
  //                  vec2.get_x(), vec2.get_y(), vec2.get_z(),
  //                  cross.get_x(), cross.get_y(), cross.get_z());
  //          m.unitarize();
  //          return m;
  //      }
  //
  //      math::Matrix
  //      _get_res_ref_frame_from_atoms(
  //              std::vector<Atom> const &atoms,
  //              ResidueTypeCOP res_type) const {
  //
  //          auto vec1 = math::Point();
  //          auto vec2 = math::Point();
  //          auto base_center = math::Point();
  //          int count = 0;
  //          for (int i = 12; i < res_type->get_num_atoms(); i++) {
  //              if (atoms[i] == nullptr) { continue; }
  //              base_center += atoms[i].get_coords();
  //              count += 1;
  //          }
  //          base_center /= float(count);
  //          auto c1p_atom = atoms[res_type->get_atom_index("C1'")];
  //
  //          if (res_type->get_short_name() == 'A' ||
  //          res_type->get_short_name() == 'G') {
  //              auto n9_atom = atoms[res_type->get_atom_index("N9")];
  //              vec1 = (n9_atom.get_coords() -
  //              c1p_atom.get_coords()).normalize(); vec2 =
  //              (n9_atom.get_coords() - base_center).normalize();
  //          } else {
  //              auto n1_atom = atoms[res_type->get_atom_index("N1")];
  //              vec1 = (n1_atom.get_coords() -
  //              c1p_atom.get_coords()).normalize(); vec2 =
  //              (n1_atom.get_coords() - base_center).normalize();
  //          }
  //
  //          auto cross = vec1.cross(vec2);
  //          auto m = math::Matrix(
  //                  vec1.get_x(), vec1.get_y(), vec1.get_z(),
  //                  vec2.get_x(), vec2.get_y(), vec2.get_z(),
  //                  cross.get_x(), cross.get_y(), cross.get_z());
  //          m.unitarize();
  //          return m;
  //      }

  // TODO Change it to private
public:
  /**
   * residue type object which explains what atoms in belong in this residue.
   */
  ResidueTypeCOP _res_type;

  /**
   * vector of the atom objects that belong to this residue
   */
  Atoms _atoms;

  /**
   * vector of bead objects for sterics
   */
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
typedef std::vector<ResidueOP> ResidueOPs;

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
///res_type->get_name() + ": IGNORING!";
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
///for residue " + spl[0] + " " + spl[1];
//              }
//          }
//      }
} // namespace structure

#endif /* defined(__RNAMake__structure_residue__) */