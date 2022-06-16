//
// Created by Joseph Yesselman on 12/17/17.
//

#ifndef RNAMAKE_NEW_SEGMENT_H
#define RNAMAKE_NEW_SEGMENT_H

#include "util/segment_type.h"
#include <primitives/pose.h>

namespace primitives {

template <typename BPtype, typename Structuretype, typename Chaintype,
          typename Restype>
class Segment : public Pose<BPtype, Structuretype, Chaintype, Restype> {
private:
  typedef Pose<BPtype, Structuretype, Chaintype, Restype> BaseClass;
  typedef std::vector<Restype> Residues;
  typedef base::VectorContainerOP<Restype> ResiduesOP;
  typedef base::VectorContainerOP<BPtype> BasepairsOP;

public:
  inline Segment(Structuretype const &structure,
                 std::vector<BPtype> const &basepairs,
                 Indexes const &end_indexes,
                 base::SimpleStringCOPs const &end_ids,
                 base::SimpleStringCOP name, util::SegmentType segment_type,
                 Index aligned_end_index, util::Uuid const &uuid)
      : BaseClass(structure, basepairs, end_indexes, end_ids, name),
        segment_type_(segment_type), aligned_end_index_(aligned_end_index) {}

  // TODO Change back to private
public:
  // let dervived classes fill in members
  Segment() : BaseClass() {}

  typedef std::vector<BPtype> Basepairs;

public:
  inline Index get_aligned_end_index() const { return aligned_end_index_; }

  util::Uuid const &get_uuid() const { return uuid_; }

  util::SegmentType get_segment_type() const { return segment_type_; }

protected:
  util::SegmentType segment_type_;
  Index aligned_end_index_;
  util::Uuid uuid_;
};

} // namespace primitives

#endif // RNAMAKE_NEW_SEGMENT_H