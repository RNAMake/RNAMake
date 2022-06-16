//
//  motif.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_motif__
#define __RNAMake__ss_motif__

#include <algorithm>
#include <stdio.h>

#include "secondary_structure/basepair.h"
#include "secondary_structure/chain.h"
#include "secondary_structure/rna_structure.h"
#include "util/motif_type.h"

namespace secondary_structure {

class Motif : public RNAStructure {
public: // RNA Structure Constructors
  Motif() : RNAStructure() {}

  Motif(StructureOP const &structure, BasepairOPs const &basepairs,
        BasepairOPs const &ends)
      : RNAStructure(structure, basepairs, ends) {}

  Motif(StructureOP const &structure, BasepairOPs const &basepairs,
        BasepairOPs const &ends, Strings const &end_ids, String const &name,
        String const &path, float score)
      : RNAStructure(structure, basepairs, ends, end_ids, name, path, score) {}

public: // Motif specific constructors
  Motif(Motif const &m) {
    structure_ = std::make_shared<Structure>(*m.structure_);
    basepairs_ = BasepairOPs(m.basepairs_.size());
    ends_ = BasepairOPs(m.ends_.size());
    int i = 0;
    for (auto const &bp : m.basepairs_) {
      auto new_bp = std::make_shared<Basepair>(
          structure_->get_residue(bp->res1()->uuid()),
          structure_->get_residue(bp->res2()->uuid()), bp->uuid());
      basepairs_[i] = new_bp;
      i++;
    }

    i = 0;
    for (auto const &end : m.ends_) {
      int pos = (int)(std::find(m.basepairs_.begin(), m.basepairs_.end(), end) -
                      m.basepairs_.begin());
      ends_[i] = basepairs_[pos];
      i++;
    }

    mtype_ = m.mtype_;
    name_ = m.name_;
    path_ = m.path_;
    end_ids_ = m.end_ids_;
  }

  Motif(String const &s) {
    auto spl = base::split_str_by_delimiter(s, "!");
    mtype_ = static_cast<util::MotifType>(std::stoi(spl[0]));
    name_ = spl[1];
    path_ = spl[2];
    structure_ = std::make_shared<Structure>(spl[3]);
    basepairs_ = BasepairOPs();
    ends_ = BasepairOPs();
    end_ids_ = base::split_str_by_delimiter(spl[6], " ");
    auto res = structure_->residues();
    for (auto const &bp_str : base::split_str_by_delimiter(spl[4], "@")) {
      auto res_is = base::split_str_by_delimiter(bp_str, " ");
      auto res1 = res[std::stoi(res_is[0])];
      auto res2 = res[std::stoi(res_is[1])];
      auto bp = std::make_shared<Basepair>(res1, res2, util::Uuid());
      basepairs_.push_back(bp);
    }
    for (auto const &end_i : base::split_str_by_delimiter(spl[5], " ")) {
      ends_.push_back(basepairs_[std::stoi(end_i)]);
    }
  }

  ~Motif() {}

public:
  String to_str() {
    auto s = String("");
    s.append(std::to_string((int)mtype_) + "!" + name_ + "!" + path_ + "!");
    s.append(structure_->to_str() + "!");
    auto res = structure_->residues();

    for (auto const &bp : basepairs_) {
      int res1_pos =
          (int)(std::find(res.begin(), res.end(), bp->res1()) - res.begin());
      int res2_pos =
          (int)(std::find(res.begin(), res.end(), bp->res2()) - res.begin());
      s.append(std::to_string(res1_pos) + " " + std::to_string(res2_pos) + "@");
    }
    s.append("!");
    for (auto const &end : ends_) {
      int bp_pos = (int)(std::find(basepairs_.begin(), basepairs_.end(), end) -
                         basepairs_.begin());
      s.append(std::to_string(bp_pos) + " ");
    }
    s.append("!");
    for (auto const &ei : end_ids_) {
      s.append(ei + " ");
    }

    return s;
  }

public: // getters
  inline util::MotifType const &mtype() { return mtype_; }

  inline util::Uuid const &id() { return id_; }

public: // setters
  inline void mtype(util::MotifType const &mtype) { mtype_ = mtype; }

  inline void id(util::Uuid const &uuid) { id_ = uuid; }

private:
  util::MotifType mtype_;
  util::Uuid id_;
};

typedef std::shared_ptr<Motif> MotifOP;
typedef std::vector<MotifOP> MotifOPs;

} // namespace secondary_structure

#endif /* defined(__RNAMake__motif__) */
