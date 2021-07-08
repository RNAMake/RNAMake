//
//  residue.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__all_atom_residue__
#define __RNAMake__all_atom_residue__

#include <stdio.h>

//RNAMake Headers
#include <util/uuid.h>
#include <util/bead.h>
#include <base/types.h>
#include <primitives/residue.h>
#include <structure/atom.h>
#include <structure/residue_type.h>
#include <structure/residue_type_set.h>
#include <doctest.h>

namespace structure {

  math::Point
  center(
          Atoms const &);

  math::Point
  center(
          AtomOPs const &);


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
 *  auto a = r->get_atom("C1'");
 *  std::cout << a->coords() << std::endl; //OUTPUT: "-23.806 -50.289  86.732"
 *
 *  auto beads = r->get_beads();
 *  std::cout << beads[0]->btype() << " " << beads[0]->center() << std::endl;
 *  //OUTPUT: "1 24.027 -48.5001111111 86.368"
 *
 *  //get PDB formatted String (can also write to file with r->to_pdb("test.pdb") )
 *  std::cout << r->to_pdb_str() << std::endl; //OUTPUT -->
 *  ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00
 *  ATOM      2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00
 *  ATOM      3 C4'  G   A 103     -24.521 -48.156  86.068  1.00  0.00
 *  ATOM      4 O4'  G   A 103     -24.861 -49.568  86.118  1.00  0.00
 *  ATOM      5 C3'  G   A 103     -23.009 -48.119  86.281  1.00  0.00
 *  ATOM      6 O3'  G   A 103     -22.548 -46.872  86.808  1.00  0.00
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
      Residue(
              char name,
              int num,
              char chain_id,
              char i_code,
              ResidueTypeCOP res_type,
              Atoms const &atoms,
              util::Uuid const &uuid) :
              primitives::Residue(name, num, chain_id, i_code, uuid),
              res_type_(res_type),
              atoms_(atoms) {
          _build_beads();
      }

      /**
       * Copy constructor
       * @param   r   Residue object copying from
       */
      Residue(
              Residue const &r) :
              primitives::Residue(r.name_, r.num_, r.chain_id_, r.i_code_, r.uuid_),
              res_type_(r.res_type_),
              atoms_(r.atoms_),
              beads_(r.beads_) {}

      /**
       * Construction from String, used in reading data from files
       * @param   s   string generated from to_str()
       * @see to_str()
       */
      Residue(
              String const &s,
              structure::ResidueTypeSet const &rts) :
              primitives::Residue() {

          auto spl = base::split_str_by_delimiter(s, ",");
          res_type_ = rts.get_residue_type(spl[0]);
          name_ = spl[1][0];
          num_ = std::stoi(spl[2]);
          chain_id_ = spl[3][0];
          i_code_ = spl[4][0];
          uuid_ = util::Uuid();
          atoms_ = Atoms();
          int i = 5;
          while (i < spl.size()) {
              if(spl[i].size() <= 1){
                  i++;
                  continue;
              }
              atoms_.push_back(Atom(spl[i]));
              i++;
          }
          _build_beads();
      }

//        Residue(
//                json::JSON &j,
//                ResidueTypeSet const &rts) :
//                primitives::Residue() {
//            res_type_ = rts.get_residue_type(j["res_type"].ToString());
//            name_ = (char) j["name"].ToInt();
//            num_ = (int) j["num"].ToInt();
//            chain_id_ = (char) j["chain_id"].ToInt();
//            i_code_ = (char) j["i_code"].ToInt();
//            uuid_ = util::Uuid();
//            atoms_ = Atoms();
//            auto &atoms_json = j["atoms"];
//            for (int i = 0; i < atoms_json.size(); i++) {
//                atoms_.push_back(Atom(atoms_json[i]));
//            }
//        }


      /**
       * Empty deconstructor
       */
      ~Residue() {}

  public: //iterator stuff
      typedef Atoms::const_iterator const_iterator;

      const_iterator begin() const noexcept { return atoms_.begin(); }

      const_iterator end() const noexcept { return atoms_.end(); }

      typedef util::Beads::const_iterator bead_const_iterator;

      bead_const_iterator bead_begin() const noexcept { return beads_.begin(); }

      bead_const_iterator bead_end() const noexcept { return beads_.end(); }

  public:

      inline
      bool
      operator==(
              Residue const &r) const {
          return is_equal(r);
      }

      inline
      bool
      operator!=(
              Residue const &r) const {
          return !is_equal(r);
      }

      friend
      std::ostream &
      operator<<(
              std::ostream &stream,
              Residue const &r) {
          stream << r.num_ << r.chain_id_;
          if (r.i_code_ != ' ') { stream << "(" << r.i_code_ << ")"; }
          return stream;
      }

      inline
      int
      connected_to(Residue const &res,
              float cutoff = 5.0) const {
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
      inline
      bool
      is_equal(
              Residue const &r,
              bool check_uuid = true) const {
          if (check_uuid && uuid_ != r.uuid_) { return false; }
          if (name_ != r.name_) { return false; }
          if (num_ != r.num_) { return false; }
          if (chain_id_ != r.chain_id_) { return false; }
          if (i_code_ != r.i_code_) { return false; }
          for (int i = 0; i < atoms_.size(); i++) {
              if (atoms_[i] != r.atoms_[i]) { return false; }
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
      inline
      Atom const &
      get_atom(
              String const &name) const {
          //TODO Figure out why it shifts by 3
          int index = res_type_->get_atom_index(name) - 3;
          if (index == -1) {
              throw ResidueException("atom name: " + name + " does not exist in this residue");
          }
          return atoms_[index];
      }

      inline
      Atom const &
      get_atom(
              Index index) const {
          return atoms_[index];
      }

      inline
      math::Point const &
      get_coords(
              String const &name) const {
          return get_atom(name).get_coords();
      }

      inline
      util::Bead const &
      get_bead(
              util::BeadType bead_type) const {
          for (auto const &b : beads_) {
              if (b.get_type() == bead_type) { return b; }
          }
          throw ResidueException(
                  "bead type: " + std::to_string((int) bead_type) + " does not exist in this residue");
      }

      // TODO Temporary get rid of these

      template<typename T>
      inline
      void
      set_chain_id(T &id) {
          chain_id_ = id;
      }

      inline
      void
      set_num(int &num) {
          num_ = num;
      }

      inline
      void
      set_uuid(util::Uuid const &id) {
          uuid_ = id;
      }

      //        inline
      //        AtomOPs const &
      //        get_atoms() const {
      //            return atoms_;
      //        }

      inline
      AtomOPs
      get_atoms() {
          AtomOPs atom_ptrs;
          for (auto const &a : atoms_) {
              AtomOP atom_ptr = std::make_shared<Atom>(a);
              atom_ptrs.push_back(atom_ptr);
          }
          return atom_ptrs;
      }

      inline
      math::Point
      get_center() const { return center(atoms_); }

      inline
      size_t
      get_num_atoms() const { return atoms_.size(); }

      inline
      String const &
      get_res_name() const { return res_type_->get_name(); }

      inline
      SetType
      get_res_set_type() const { return res_type_->get_set_type(); }


      /**
       * stringifes residue object
       * @code
       *  #include "instances/structure_instances.hpp"
       *  auto r = instances::residue();
       *  std::cout << r->to_str() << std::endl;
       *  //OUTPUT
       *  GUA,G,103,A,,N,N,N,O5' -26.469 -47.756 84.669,C5' -25.05 -47.579 84.775,C4' -24.521 -48.156
       *  86.068,O4' -24.861 -49.568 86.118,C3' -23.009 -48.119 86.281,O3' -22.548 -46.872 86.808,C1'
       *  -23.806 -50.289 86.732,C2' -22.812 -49.259 87.269,O2' -23.167 -48.903 88.592,N1 -19.538 -52.485
       *  85.025,C2 -19.717 -51.643 86.097,N2 -18.624 -51.354 86.809,N3 -20.884 -51.124 86.445,C4 -21.881
       *  -51.521 85.623,C5 -21.811 -52.356 84.527,C6 -20.546 -52.91 84.164,O6 -20.273 -53.677 83.228,N7
       *  -23.063 -52.513 83.947,C8 -23.858 -51.786 84.686,N9 -23.21 -51.159 85.722,
       * @endcode
       */
      String
      get_str() const;

//      json::JSON
//      get_json() const;

      /**
       * wrapper for to_pdb_str(int &, int, String const &) when one does not care about
       * renumbering atoms and residue
       **/
      inline
      String
      get_pdb_str(
              int acount) const {
          auto num = num_;
          auto chain_id = chain_id_;
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
       *  //get PDB formatted String (can also write to file with r->to_pdb("test.pdb") )
       *  std::cout << r->to_pdb_str() << std::endl; //OUTPUT -->
       *  ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00
       *  ATOM      2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00
       *  ATOM      3 C4'  G   A 103     -24.521 -48.156  86.068  1.00  0.00
       *  ATOM      4 O4'  G   A 103     -24.861 -49.568  86.118  1.00  0.00
       *  ATOM      5 C3'  G   A 103     -23.009 -48.119  86.281  1.00  0.00
       *  ATOM      6 O3'  G   A 103     -22.548 -46.872  86.808  1.00  0.00
       *  .
       *  .
       *  .
       * @endcode
       */
      String
      get_pdb_str(
              int &,
              int,
              char) const;

      String
      get_bead_pdb_str(
              int &,
              int,
              char) const;

      /**
       * writes a PDB string formmated verision of this Residue object to file
       * @param  filename of output PDB file
       */
      void
      write_pdb(String const);

  public: // non const

      inline
      void
      move(
              math::Point const &p) {
          for (auto &a : atoms_) { a.move(p); }
          for (auto &b : beads_) { b.move(p); }
      }

      inline
      void
      transform(
              math::Matrix const &r,
              math::Vector const &t,
              math::Point &dummy) {
          for (auto &a : atoms_) { a.transform(r, t, dummy); }
          for (auto &b : beads_) { b.transform(r, t, dummy); }
      }

      inline
      void
      transform(
              math::Matrix const &r,
              math::Vector const &t) {
          auto dummy = math::Point();
          transform(r, t, dummy);
      }

      inline
      void
      remove_beads() { beads_ = util::Beads(); }

      inline
      void
      build_beads() { _build_beads(); }

      inline
      void
      new_uuid() { uuid_ = util::Uuid(); }


  public: // getters

      /**
       * getter for the one letter residue type
       */
      /*inline
      String
      short_name() const { return res_type_->short_name(); }*/

  private:
      void
      _build_beads();

      void
      _build_beads_RNA();


      //TODO Change it to private
  public:
      /**
       * residue type object which explains what atoms in belong in this residue.
       */
      ResidueTypeCOP res_type_;

      /**
       * vector of the atom objects that belong to this residue
       */
      Atoms atoms_;


      /**
       * vector of bead objects for sterics
       */
      util::Beads beads_;

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

  inline
  bool
  residue_steric_clash_RNA(
          Residue const &r1,
          Residue const &r2) {
      for (auto it1 = r1.bead_begin(); it1 != r1.bead_end(); it1++) {
          if (it1->get_type() == util::BeadType::PHOS) { continue; }
          for (auto it2 = r2.bead_begin(); it2 != r2.bead_end(); it2++) {
              if (it2->get_type() == util::BeadType::PHOS) { continue; }
              if (it1->distance(*it2) < 2.5) { return true; }
          }
      }

      /*std::for_each(r1.bead_begin(), r1.bead_end(), [](util::Bead const & b1) {
          if(b1.get_type() != util::BeadType::PHOS) {
              std::for_each(r2.bead_begin(), r2.bead_end(), [](util::Bead const & b2) {
                  if (b2.get_type() != util::BeadType::PHOS) {
                      if (b1.distance(b2) < 2.5) { return true; }
                  }
              });
          }
      });*/

      return false;
  }

  inline
  bool
  residue_steric_clash(
          Residue const &r1,
          Residue const &r2) {
      for (auto it1 = r1.bead_begin(); it1 != r1.bead_end(); it1++) {
          for (auto it2 = r2.bead_begin(); it2 != r2.bead_end(); it2++) {
              if (it1->distance(*it2) < 2.5) { return true; }
          }
      }

      return false;
  }
}

#endif /* defined(__RNAMake__all_atom_residue__) */