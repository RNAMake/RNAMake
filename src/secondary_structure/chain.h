
//
//  chain.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_chain__
#define __RNAMake__ss_chain__

#include <stdio.h>

#include <primitives/chain.h>
#include <secondary_structure/residue.h>

namespace secondary_structure {

class Chain : public primitives::Chain<Residue> {
public:
  typedef primitives::Chain<Residue> ParentClass;

public:
  inline Chain(Residues const &residues) : ParentClass(residues) {}

  inline Chain(Chain const &c) : ParentClass(c) {}

  inline Chain(String const &s) : ParentClass(s) {}

  virtual ~Chain() {}

public:
  inline bool operator==(Chain const &c) const {
    if (residues_.size() != c.residues_.size()) {
      return false;
    }

    for (int i = 0; i < c.get_length(); i++) {
      if (residues_[i] != c.residues_[i]) {
        return false;
      }
    }
    return true;
  }

  inline bool operator!=(Chain const &c) const { return !(*this == c); }

public:
  bool is_equal(Chain const &c, bool check_uuid = true) const {

    if (residues_.size() != c.residues_.size()) {
      return false;
    }

    for (int i = 0; i < c.get_length(); i++) {
      if (residues_[i].is_equal(c.residues_[i]), check_uuid) {
        return false;
      }
    }
    return true;
  }

public:
  inline String get_dot_bracket() const {
    auto db = String("");
    for (auto const &r : residues_) {
      db += r.get_dot_bracket();
    }
    return db;
  }

  inline String get_sequence() const {
    auto seq = String("");
    for (auto const &r : residues_) {
      seq += r.get_name();
    }
    return seq;
  }

  inline String get_str() const {
    auto s = String("");
    for (auto const &r : residues_) {
      s += r.get_str() + ";";
    }
    return s;
  }

public:
  inline void set_residue_identity(Index residue_index, char name) {
    residues_[residue_index].set_name(name);
  }
};

typedef std::shared_ptr<Chain> ChainOP;
typedef std::vector<ChainOP> ChainOPs;

} // namespace secondary_structure

#endif /* defined(__RNAMake__chain__) */