//
// Created by Joe Yesselman on 6/28/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_STATE_SEGMENT_HPP_
#define RNAMAKE_SRC_STRUCTURE_STATE_SEGMENT_HPP_

#include <structure/base/segment.hpp>
#include <structure/state/residue.hpp>
#include <structure/state/basepair.hpp>

namespace structure::state {

typedef structure::base::Chain<Residue> Chain;
typedef structure::base::Structure<Chain, Residue> Structure;
typedef structure::base::Pose<Basepair, Structure, Chain, Residue> Pose;
typedef structure::base::Segment<Basepair, Structure, Chain, Residue> Segment;

}

#endif // RNAMAKE_SRC_STRUCTURE_STATE_SEGMENT_HPP_
