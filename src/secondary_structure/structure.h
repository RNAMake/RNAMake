
//
//  structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sec_structure__
#define __RNAMake__sec_structure__

#include <stdio.h>

// RNAMake Headers
#include <base/types.hpp>
#include <primitives/structure.h>
#include <secondary_structure/chain.h>

namespace secondary_structure {

class Structure : public primitives::Structure<Chain, Residue> {
public:
  typedef primitives::Structure<Chain, Residue> ParentClass;

public:
  inline Structure(Residues const &residues, Cutpoints const &cut_points)
      : ParentClass(residues, cut_points) {}
  inline Structure(Structure const &structure) : ParentClass(structure) {}

  inline Structure(String const &s) : ParentClass(s) {}

  ~Structure() {}

public:
  bool operator==(Structure const &s) const {
    if (_residues.size() != s._residues.size()) {
      return false;
    }
    if (_cut_points.size() != s._cut_points.size()) {
      return false;
    }

    for (int i = 0; i < _residues.size(); i++) {
      if (_residues[i] != s._residues[i]) {
        return false;
      }
    }

    for (int i = 0; i < _cut_points.size(); i++) {
      if (_cut_points[i] != s._cut_points[i]) {
        return false;
      }
    }

    return true;
  }

  bool operator!=(Structure const &s) const { return !(*this == s); }

public:
  bool is_equal(Structure const &s, bool check_uuid = true) const {

    if (_residues.size() != s._residues.size()) {
      return false;
    }
    if (_cut_points.size() != s._cut_points.size()) {
      return false;
    }

    for (int i = 0; i < _residues.size(); i++) {
      if (_residues[i].is_equal(_residues[i]), check_uuid) {
        return false;
      }
    }

    for (int i = 0; i < _cut_points.size(); i++) {
      if (_cut_points[i] != s._cut_points[i]) {
        return false;
      }
    }
    return true;
  }

public: // getters
  String get_str() {
    auto s = String();
    for (auto const &r : _residues) {
      s += r.get_str() + ";";
    }
    int i = 0;
    for (auto const &c : _cut_points) {
      s += std::to_string(c);
      if (i != _cut_points.size()) {
        s += " ";
      }
    }
    s += ";";
    return s;
  }

  String get_dot_bracket() const {
    auto i = -1;
    auto seq = String("");
    auto pos = 0;
    for (auto const &r : _residues) {
      i++;
      if (_cut_points[pos] == i) {
        seq += "&";
        pos++;
      }
      seq += r.get_dot_bracket();
    }
    return seq;
  }

public: // setters
  void set_sequence(String const &sequence) {
    auto spl = base::split_str_by_delimiter(sequence, "&");
    auto new_sequence = base::join_by_delimiter(spl, "");
    expects<StructureException>(new_sequence.length() == _residues.size(),
                                "invalid length of new sequence, has a "
                                "different number then number of residues");

    auto i = 0;
    for (auto &r : _residues) {
      r.set_name(new_sequence[i]);
      i++;
    }
  }

  inline void set_residue_identity(Index residue_index, char name) {
    _residues[residue_index].set_name(name);
  }
};

typedef std::shared_ptr<Structure> StructureOP;
typedef std::shared_ptr<Structure> StructureUP;

StructureOP get_structure_from_secondary_structure(String const &,
                                                   String const &);

} // namespace secondary_structure

#endif /* defined(__RNAMake__structure__) */