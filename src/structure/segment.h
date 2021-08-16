
//
// Created by Joseph Yesselman on 12/20/17.
//

#ifndef RNAMAKE_NEW_STRUCTURE_SEGMENT_H
#define RNAMAKE_NEW_STRUCTURE_SEGMENT_H

#include <util/segment_type.h>
#include <primitives/segment.h>
#include <secondary_structure/segment.h>
#include <structure/structure.h>
#include <structure/basepair.h>

// forward declaration to allow SegmentSqliteLibrary to call new_uuids()
namespace resources {
  class SegmentSqliteLibrary;
}

    namespace structure {

  class Segment : public primitives::Segment<Basepair, Structure, Chain, Residue>  {
  public:
      typedef primitives::Segment<Basepair, Structure, Chain, Residue> BaseClass;

  public:
      friend class SegmentFactory;
      friend class resources::SegmentSqliteLibrary;

  public:
      Segment(
              Structure const & structure,
              Basepairs const & basepairs,
              Indexes const & end_indexes,
              base::SimpleStringCOPs const & end_ids,
              base::SimpleStringCOP name,
              Structure const & proteins,
              Structure const & small_molecules,
              base::SimpleStringCOP dot_bracket,
              util::SegmentType segment_type,
              Index aligned_end_index,
              util::Uuid const & uuid ):
              BaseClass(structure, basepairs, end_indexes, end_ids, name, segment_type, aligned_end_index, uuid),
              proteins_(proteins),
              small_molecules_(small_molecules),
              dot_bracket_(dot_bracket) {}

      Segment(
              Segment const & seg):
              BaseClass(seg.structure_, seg.basepairs_, seg.end_indexes_, seg.end_ids_, seg.name_,
                        seg.segment_type_, seg.aligned_end_index_, seg.uuid_),
              proteins_(seg.proteins_),
              small_molecules_(seg.small_molecules_),
              dot_bracket_(seg.dot_bracket_) {}


  public:

      const_iterator protein_begin() const { return proteins_.begin(); }
      const_iterator protein_end()   const { return proteins_.end(); }

      inline
      bool
      is_protein_residue_start_of_chain(
              Residue const & r) const {
          return proteins_.is_residue_start_of_chain(r);
      }

      const_iterator small_molecules_begin() const { return small_molecules_.begin(); }
      const_iterator small_molecules_end()   const { return small_molecules_.end(); }

  public:
      bool
      is_equal(
              Segment const & s,
              bool check_uuid = true) const {
          if(segment_type_ != s.segment_type_) { return false; }
          if(aligned_end_index_ != s.aligned_end_index_) { return false; }
          if(basepairs_.size() != s.basepairs_.size()) { return false; }
          if(end_indexes_.size() != s.end_indexes_.size()) { return false; }
          if(*name_ != *s.name_) { return false; }
          if(*dot_bracket_ != *s.dot_bracket_) { return false; }
          for(int i = 0; i < basepairs_.size(); i++) {
              if(!basepairs_[i].is_equal(s.basepairs_[i], check_uuid)) { return false; }
          }
          for(int i = 0; i < end_indexes_.size(); i++) {
              if(end_indexes_[i] != s.end_indexes_[i]) { return false; }
          }
          if(! structure_.is_equal(s.structure_, check_uuid)) { return false; }
          if(! proteins_.is_equal(s.proteins_, check_uuid)) { return false; }
          if(! small_molecules_.is_equal(s.small_molecules_, check_uuid)) { return false; }
          return true;
      }



  public: // non const methods
      void
      move(
              math::Point const & p) {
          structure_.move(p);
          proteins_.move(p);
          small_molecules_.move(p);
          for(auto & bp : basepairs_) { bp.move(p); }

      }

      void
      transform(
              math::Matrix const & r,
              math::Vector const & t,
              math::Point & dummy) {
          structure_.transform(r, t, dummy);
          proteins_.transform(r, t, dummy);
          small_molecules_.transform(r, t, dummy);
          for(auto & bp : basepairs_) { bp.transform(r, t, dummy); }
      }

      inline
      void
      transform(
              math::Matrix const & r,
              math::Vector const & t) {
          auto dummy = math::Point();
          transform(r, t, dummy);
      }

  public:
      bool
      steric_clash(
              Segment const & s) {

          for(auto const & r1 : *this) {
              // RNA/RNA clashes
              for(auto const & r2 : s) {
                  if(residue_steric_clash_RNA(r1, r2)) { return true; }
              }
              // RNA/protein clashes
              for(auto const & r2 : s.proteins_) {
                  if(residue_steric_clash_RNA(r1, r2)) { return true; }
              }
              // RNA/small molecule clashes
              for(auto const & r2 : s.small_molecules_) {
                  if(residue_steric_clash_RNA(r1, r2)) { return true; }
              }
          }

          for(auto const & r1 : proteins_) {
              // protein/RNA clashes
              for(auto const & r2 : s) {
                  if(residue_steric_clash_RNA(r1, r2)) { return true; }
              }
              // protein/protein clashes
              for(auto const & r2 : s.proteins_) {
                  if(residue_steric_clash(r1, r2)) { return true; }
              }
              // protein/small molecule clashes
              for(auto const & r2 : s.small_molecules_) {
                  if(residue_steric_clash(r1, r2)) { return true; }
              }
          }

          for(auto const & r1 : small_molecules_) {
              // small_molecule/RNA clashes
              for(auto const & r2 : s) {
                  if(residue_steric_clash_RNA(r1, r2)) { return true; }
              }
              // small_molecule/protein clashes
              for(auto const & r2 : s.proteins_) {
                  if(residue_steric_clash(r1, r2)) { return true; }
              }
              // small_molecule/small molecule clashes
              for(auto const & r2 : s.small_molecules_) {
                  if(residue_steric_clash(r1, r2)) { return true; }
              }
          }
          return false;
      }

      int
      get_num_steric_clashes(
              Segment const & s) {
          int steric_clash_count = 0;

          for(auto const & r1 : *this) {
              // RNA/RNA clashes
              for(auto const & r2 : s) {
                  if(residue_steric_clash_RNA(r1, r2)) { steric_clash_count += 1; }
              }
              // RNA/protein clashes
              for(auto const & r2 : s.proteins_) {
                  if(residue_steric_clash_RNA(r1, r2)) { steric_clash_count += 1;  }
              }
              // RNA/small molecule clashes
              for(auto const & r2 : s.small_molecules_) {
                  if(residue_steric_clash_RNA(r1, r2)) { steric_clash_count += 1;  }
              }
          }

          for(auto const & r1 : proteins_) {
              // protein/RNA clashes
              for(auto const & r2 : s) {
                  if(residue_steric_clash_RNA(r1, r2)) { steric_clash_count += 1;  }
              }
              // protein/protein clashes
              for(auto const & r2 : s.proteins_) {
                  if(residue_steric_clash(r1, r2)) { steric_clash_count += 1;  }
              }
              // protein/small molecule clashes
              for(auto const & r2 : s.small_molecules_) {
                  if(residue_steric_clash(r1, r2)) { steric_clash_count += 1;  }
              }
          }

          for(auto const & r1 : small_molecules_) {
              // small_molecule/RNA clashes
              for(auto const & r2 : s) {
                  if(residue_steric_clash_RNA(r1, r2)) { steric_clash_count += 1;  }
              }
              // small_molecule/protein clashes
              for(auto const & r2 : s.proteins_) {
                  if(residue_steric_clash(r1, r2)) { steric_clash_count += 1;  }
              }
              // small_molecule/small molecule clashes
              for(auto const & r2 : s.small_molecules_) {
                  if(residue_steric_clash(r1, r2)) { steric_clash_count += 1;  }
              }
          }

          return steric_clash_count;
      }


      secondary_structure::SegmentOP
      get_secondary_structure() const {
          auto dot_bracket_spl = base::split_str_by_delimiter(dot_bracket_->get_str(), "&");
          auto dot_bracket_str = base::join_by_delimiter(dot_bracket_spl, "");

          auto i = 0 ;
          auto residues = secondary_structure::Residues();
          auto cutpoints = structure_.get_cutpoints();
          for(auto const & r : structure_) {
              auto ss_r = secondary_structure::Residue(r.get_name(), dot_bracket_str[i], r.get_num(),
                                                       r.get_chain_id(), r.get_i_code(), r.get_uuid());
              residues.push_back(ss_r);
              i += 1;
          }

          auto ss_structure = secondary_structure::Structure(residues, cutpoints);
          auto ss_bps = secondary_structure::Basepairs();
          for(auto const & bp : basepairs_) {
              auto ss_bp = secondary_structure::Basepair(
                      bp.get_res1_uuid(), bp.get_res2_uuid(), bp.get_uuid(), bp.get_bp_type(), bp.get_name());

              ss_bps.push_back(ss_bp);
          }

          return std::make_shared<secondary_structure::Segment>(
                  ss_structure, ss_bps, end_indexes_, end_ids_, name_, segment_type_, aligned_end_index_, uuid_);
      }


  public: // pdb functions
      String
      get_pdb_str(
              int &,
              int &,
              char &);

      inline
      String
      get_pdb_str(
              int acount = 0) {
          auto num = structure_.get_residue(0).get_num();
          auto chain_id = structure_.get_residue(0).get_chain_id();
          return get_pdb_str(acount, num, chain_id);
      }

      void
      write_pdb(
              String const &) const;

      void
      write_steric_beads_to_pdb(
              String const &);

  public:
      inline
      String
      get_dot_bracket() { return dot_bracket_->get_str(); }

      inline
      Basepair const &
      get_aligned_end() { return basepairs_[end_indexes_[aligned_end_index_]]; }

  protected:
      void
      new_uuids() {
          uuid_ = util::Uuid();
          auto uuid_map = std::map<util::Uuid, int>();
          int i = 0;
          for(auto const & r : structure_) {
              uuid_map[r.get_uuid()] = i;
              i++;
          }
          structure_.new_uuids();

          for(auto & bp : basepairs_) {
              auto r1_pos = uuid_map[bp.get_res1_uuid()];
              auto r2_pos = uuid_map[bp.get_res2_uuid()];
              auto & r1 = structure_.get_residue(r1_pos);
              auto & r2 = structure_.get_residue(r2_pos);
              bp.new_uuids(r1.get_uuid(), r2.get_uuid());
          }

          proteins_.new_uuids();
          small_molecules_.new_uuids();
      }

  protected:
      Structure proteins_;
      Structure small_molecules_;
      base::SimpleStringCOP dot_bracket_;


  };

  typedef std::shared_ptr<Segment> SegmentOP;
  typedef std::vector<SegmentOP>   SegmentOPs;


}


#endif //RNAMAKE_STRUCTURE_SEGMENT_H