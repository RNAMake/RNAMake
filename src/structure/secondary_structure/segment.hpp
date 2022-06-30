//
// Created by Joe Yesselman on 6/27/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_SECONDARY_STRUCTURE_SEGMENT_HPP_
#define RNAMAKE_SRC_STRUCTURE_SECONDARY_STRUCTURE_SEGMENT_HPP_

#include <structure/secondary_structure/residue.h>
#include <structure/secondary_structure/basepair.h>

namespace structure::secondary_structure {
typedef structure::Chain<Residue> Chain;
typedef structure::Structure<Chain, Residue> Structure;
typedef structure::Pose<Basepair, Structure, Chain, Residue> Pose;
typedef structure::Segment<Basepair, Structure, Chain, Residue> Segment;

}

#endif // RNAMAKE_SRC_STRUCTURE_SECONDARY_STRUCTURE_SEGMENT_HPP_
