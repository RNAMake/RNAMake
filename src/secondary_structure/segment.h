//
// Created by Joseph Yesselman on 11/3/17.
//
#ifndef RNAMAKE_NEW_RNA_SEGMENT_H
#define RNAMAKE_NEW_RNA_SEGMENT_H

#include <util/segment_type.h>
#include <primitives/segment.h>
#include <secondary_structure/pose.h>

namespace secondary_structure {

  class Segment : public primitives::Segment<Basepair, Structure, Chain, Residue>  {
  public:
      typedef primitives::Segment<Basepair, Structure, Chain, Residue> BaseClass;

  public:
      Segment(
              Structure const & structure,
              Basepairs const & basepairs,
              Indexes const & end_indexes,
              base::SimpleStringCOPs const & end_ids,
              base::SimpleStringCOP name,
              util::SegmentType segment_type,
              Index aligned_end_index,
              util::Uuid const & uuid):
              BaseClass(structure, basepairs, end_indexes, end_ids, name, segment_type, aligned_end_index, uuid),
              sequence_update_(false) {}

  public:

      inline
      String
      get_dot_bracket() { return structure_.get_dot_bracket(); }

      inline
      void
      set_sequence(
              String const & sequence) {
          structure_.set_sequence(sequence);
          sequence_update_ = true;
      }

      inline
      void
      set_residue_identity(
              Index residue_index,
              char name) {
          structure_.set_residue_identity(residue_index, name);
          sequence_update_ = true;
      }

  public: //overrided

      base::SimpleStringCOP
      get_end_id(
              Index index) const{

          if(sequence_update_) {
              _update_end_ids();
              sequence_update_ = false;
          }

          return primitives::Segment<Basepair, Structure, Chain, Residue>::get_end_id(index);
      }


  private:
      void
      _update_end_ids() const {
          int i = -1;
          for(auto const & ei : this->end_indexes_) {
              i++;
              auto end_id = generate_end_id(this->structure_, this->basepairs_, this->basepairs_[ei]);
              this->end_ids_[i] = std::make_shared<base::SimpleString const>(end_id);
          }

      }


  private:
      mutable bool sequence_update_;


  };

  typedef std::shared_ptr<Segment> SegmentOP;

}



#endif //RNAMAKE_NEW_RNA_SEGMENT_H