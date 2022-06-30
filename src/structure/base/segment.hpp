//
// Created by Joe Yesselman on 6/30/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_BASE_SEGMENT_HPP_
#define RNAMAKE_SRC_STRUCTURE_BASE_SEGMENT_HPP_

#include <structure/base/pose.hpp>

namespace structure::base {

template <typename Basepair, typename Structure, typename Chain,
          typename Residue>
class Segment : public Pose<Basepair, Structure, Chain, Residue> {
private:
  typedef Pose<Basepair, Structure, Chain, Residue> BaseClass;
  typedef std::vector<Residue> Residues;
  typedef std::vector<Basepair> Basepairs;

public:
  inline Segment(Structure &structure, Structure &proteins,
                 Structure &small_molecules, Basepairs &basepairs,
                 Indexes &end_indexes, Strings &end_ids, String &name,
                 String &dot_bracket, util::MotifType segment_type,
                 Index aligned_end_index, util::Uuid const &uuid)
      : BaseClass(structure, proteins, small_molecules, basepairs, end_indexes,
                  end_ids, name, dot_bracket),
        _segment_type(segment_type), _aligned_end_index(aligned_end_index),
        _uuid(uuid) {}

public: // trival getters ////////////////////////////////////////////////////
  [[nodiscard]] inline Index get_aligned_end_index() const {
    return _aligned_end_index;
  }

  [[nodiscard]] inline util::Uuid const &get_uuid() const { return _uuid; }

  [[nodiscard]] inline util::MotifType get_segment_type() const {
    return _segment_type;
  }

  [[nodiscard]] inline const Basepair &get_aligned_end() const {
    return this->_basepairs[this->_end_indexes[_aligned_end_index]];
  }

private:
  util::MotifType _segment_type;
  Index _aligned_end_index;
  util::Uuid _uuid;
};

}

#endif // RNAMAKE_SRC_STRUCTURE_BASE_SEGMENT_HPP_
