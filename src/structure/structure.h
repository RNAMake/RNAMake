//
// Created by Joseph Yesselman on 11/28/17.
//

#ifndef RNAMAKE_NEW_STRUCTURE_H
#define RNAMAKE_NEW_STRUCTURE_H

#include <base/vector_container.h>
#include <util/x3dna.h>
#include <primitives/structure.h>
#include <structure/residue.h>
#include <structure/chain.h>
#include <structure/basepair.h>

namespace structure {

  class Structure : public primitives::Structure<Chain, Residue> {
  public:
      typedef primitives::Structure<Chain, Residue> ParentClass;

  public:

      inline
      Structure() : ParentClass() {}

      inline
      Structure(
              Residues const & residues,
              Cutpoints const & cut_points): ParentClass(residues, cut_points) {}
      inline
      Structure(
              Structure const & structure): ParentClass(structure._residues, structure._cut_points) {}
      inline
      Structure(
              String const & s,
              ResidueTypeSet const & rts) : ParentClass() {

          auto spl = base::split_str_by_delimiter(s, ";");
          for(Index i = 0; i < spl.size()-1; i++) {
              _residues.push_back(Residue(spl[i], rts));
          }
//          auto cut_point_spl = base::split_str_by_delimiter(spl.back(), " ");
//          for(auto const & cut_point_s : cut_point_spl) {
//              _cut_points.push_back(std::stoi(cut_point_s));
//          }
      }


      ~Structure() {}

  public: //operators
      inline
      bool
      operator == (
              Structure const & s) const {
          return is_equal(s);
      }

      inline
      bool
      operator !=(
              Structure const & s) const {
          return !(is_equal(s));
      }

  public:
      bool
      is_equal(
              Structure const & s,
              bool check_uuid = true) const {

          if(_residues.size() != s._residues.size() ) { return false; }
          if(_cut_points.size() != s._cut_points.size() ) { return false; }

          for(int i = 0; i < _residues.size(); i++) {
              if(!(_residues[i].is_equal(s._residues[i], check_uuid))) { return false; }
          }

          for(int i = 0; i < _cut_points.size(); i++) {
              if(_cut_points[i] != s._cut_points[i]) { return false; }
          }
          return true;

      }

  public: // non const
      void
      move(
              math::Point const & p) {
          for(auto & r : _residues) { r.move(p); }
      }

      void
      transform(
              math::Matrix const & r,
              math::Vector const & t,
              math::Point & dummy) {

          for(auto & res : _residues) { res.transform(r, t, dummy); }
      }

      inline
      void
      transform(
              math::Matrix const & r,
              math::Vector const & t) {
          auto dummy = math::Point();
          transform(r, t, dummy);
      }


      inline
      void
      remove_residue_beads(
              util::Uuid const & r_uuid) {

          auto & r = get_residue(r_uuid);
          auto i =  get_res_index(r);
          _residues[i].remove_beads();
      }

      void
      new_uuids() {
          for(auto & r : _residues) { r.new_uuid(); }
      }

  public: //getters

      String
      get_str() const{
          auto s = String();
          for(auto const & r : _residues) {
              s += r.get_str() + ";";
          }
          int i = 0;
          for(auto const & c : _cut_points) {
              s += std::to_string(c);
              if(i != _cut_points.size()) { s += " "; }
          }
          s += ";";
          return s;
      }


      String
      get_pdb_str(
              int &,
              int &,
              char &) const;

      inline
      String
      get_pdb_str(
              int acount = 0) const {
          auto num = _residues[0].get_num();
          auto chain_id = _residues[0].get_chain_id();
          return get_pdb_str(acount, num, chain_id);
      }

      void
      write_pdb(
              String const &) const;

      void
      write_steric_beads_to_pdb(
              String const &);


  };

  typedef std::shared_ptr<Structure> StructureOP;

  int
  are_residues_connected_RNA(
          Residue const &,
          Residue const &);

  int
  are_residues_connected_protein(
          Residue const &,
          Residue const &);

  int
  are_residues_connected(
          Residue const &,
          Residue const &);

  ResidueOP
  _get_first_residues_in_chain(
          ResidueOPs const &,
          std::map<ResidueOP, int> const &);

  ResidueOP
  _get_next_residue(
          ResidueOPs const &,
          std::map<ResidueOP, int> const &);

  StructureOP
  get_structure_from_residues(
          ResidueOPs const &);

  StructureOP
  get_structure_from_pdb(
          String const &,
          ResidueTypeSet const &,
          SetType);

  base::VectorContainerOP<Basepair>
  get_basepairs_from_x3dna(
          util::X3dna::X3Basepairs const &,
          Structure const &);


}


#endif //RNAMAKE_NEW_STRUCTURE_H