//
// Created by Joe Yesselman on 6/27/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_SECONDARY_STRUCTURE_SEGMENT_HPP_
#define RNAMAKE_SRC_STRUCTURE_SECONDARY_STRUCTURE_SEGMENT_HPP_

#include <structure/base/segment.hpp>
#include <structure/secondary_structure/residue.h>
#include <structure/secondary_structure/basepair.h>

namespace structure::secondary_structure {
  typedef structure::base::Chain<Residue> Chain;
  typedef structure::base::Structure<Chain, Residue> Structure;
  typedef structure::base::Pose<Basepair, Structure, Chain, Residue> Pose;
  typedef structure::base::Segment<Basepair, Structure, Chain, Residue> Segment;
}

#endif // RNAMAKE_SRC_STRUCTURE_SECONDARY_STRUCTURE_SEGMENT_HPP_
