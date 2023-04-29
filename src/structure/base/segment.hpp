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

  String to_str() const {
    // Make sure all of the methods are const methods
    String s = String("&");
    // Skipping path
    s += this->_name; // e.g., "HELIX.IDEAL.2"
    s += "&";
    // Need to add ???? plus &... actually, I don't think it's used
    s += "&";
    s += std::to_string(get_aligned_end_index());
    s += "&";
    s += std::to_string((int)get_segment_type());
    s += "&";
    // Structure to string:
    s += this->_structure.get_str();
    // More detailed basepair to string (?):
    s += "&";
    s += basepair_strings();
    // Basepairs to string:
    s += "&";
    // Ends to string:
    for (Index end_index : this->_end_indexes) {
      s += std::to_string(end_index) + " ";
    }
    s += "&";
    // End IDs to string
    for (auto end_id : this->_end_ids) {
      s += end_id;
      s += " ";
    }
    s += "&";
    s += secondary_structure_to_str();
    return s;
  }

  String secondary_structure_to_str() const {
    String s = String("");
    s += "99!assembled!assembled!"; // Hardcoded to match old code, isn't required
    int index = 0;
    String dot_bracket = this->get_dot_bracket();
    for (auto residue : this->_structure.get_residues()) {
      auto bracket = dot_bracket[index];
      if (bracket == '&') {
        s += "|";
        index++;
        bracket = dot_bracket[index];
      }
      s += residue.get_name();
      s += ",";
      s += bracket;
      s += ",";
      s += std::to_string(residue.get_num());
      s += ",";
      s += residue.get_chain_id();
      s += ",;";
      index++;
    }
    s += "|!";
    s += bp_to_str();
    s += "!";
    // Ends to string:
    for (Index end_index : this->_end_indexes) {
      s += std::to_string(end_index);
      s += " ";
    }
    s += "!";
    // End IDs to string
    for (auto end_id : this->_end_ids) {
      s += end_id;
      s += " ";
    }
    s += "&&";
    return s;
  }

  String basepair_strings() const {
    String s = String("");
    for (auto bp : this->_basepairs) {
      s += bp.get_name() + "," + bp.get_str() + "@";
    }
    return s;
  }

  String bp_to_str() const {
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

  inline bool operator==(const Segment &s) {
    return true;
  }

private:
  util::MotifType _segment_type;
  Index _aligned_end_index;
  util::Uuid _uuid;
};

}

#endif // RNAMAKE_SRC_STRUCTURE_BASE_SEGMENT_HPP_
