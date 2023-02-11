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

  /// @brief - gets the UUID of the segment in the structure
  [[nodiscard]] inline util::Uuid const &get_uuid() const { return _uuid; }

  /// @brief - gets the type of the segment in the structure
  [[nodiscard]] inline util::MotifType get_segment_type() const {
    return _segment_type;
  }

  /// @brief - ???
  [[nodiscard]] inline const Basepair &get_aligned_end() const {
    return this->_basepairs[this->_end_indexes[_aligned_end_index]];
  }

  String ss_to_str() {
    String s = String("");
    s += std::to_string((int)_segment_type); // e.g., "99"
    s += "!";
    s += this->_name; // e.g., "assembled"
    s += "!";
    // Structure to string:
    s += this->_structure.get_str();
    s += "|";
    // Basepairs to string:
    s += bp_to_str();
    s += "!";
    // Ends to string:
    //
    s += "!";
    // End IDs to string
    //
    s += "!";
    return s;
  }

  String bp_to_str() {
    String s = String("");
    for (auto bp : this->_basepairs) {
      auto res1_uuid = bp.get_res1_uuid();
      auto res2_uuid = bp.get_res2_uuid();
      int res1_pos, res2_pos;
      int cursor = 0;
      for (auto residue : this->_structure.get_residues()) {
        if (residue.get_uuid() == res1_uuid) { res1_pos = cursor; }
        if (residue.get_uuid() == res2_uuid) { res2_pos = cursor; }
        cursor++;
      }
      s += std::to_string(res1_pos) + " " + std::to_string(res2_pos) + "@";
    }
    return s;
  }

private:
  util::MotifType _segment_type;
  Index _aligned_end_index;
  util::Uuid _uuid;
};

}

#endif // RNAMAKE_SRC_STRUCTURE_BASE_SEGMENT_HPP_
