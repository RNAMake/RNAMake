//
//  pose.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__pose__
#define __RNAMake__pose__

#include <stdio.h>

#include <map>
#include <memory>

// RNAMake Headers
#include "motif/motif.h"
#include "util/motif_type.h"

namespace motif {

class Pose : public Motif {
 public:
  Pose()
      : designable_(std::map<util::Uuid, int, util::UuidCompare>()),
        motifs_(std::map<util::MotifType, motif::MotifOPs>()) {}

  Pose(MotifOP const &);

  Pose(structure::StructureOP const &, structure::BasepairOPs const &);

  ~Pose() {}

 public:
  motif::MotifOPs const &motifs(util::MotifType const &);

 public:
  void set_bp_designable(structure::BasepairOP const &);

  String designable_sequence();

 public:
  inline void designable(
      std::map<util::Uuid, int, util::UuidCompare> const &ndesignable) {
    designable_ = ndesignable;
  }

  inline void set_motifs(
      std::map<util::MotifType, motif::MotifOPs> const &motifs) {
    motifs_ = motifs;
  }

  /*inline
  void
  set_ss_motifs(std::map<String, secondary_structure::motif::MotifOPs> &
  ss_motifs) { secondary_structure_->set_motifs(ss_motifs);
  }*/

 private:
  std::map<util::Uuid, int, util::UuidCompare> designable_;
  std::map<util::MotifType, motif::MotifOPs> motifs_;
};

typedef std::shared_ptr<Pose> PoseOP;

}  // namespace motif

#endif /* defined(__RNAMake__pose__) */
