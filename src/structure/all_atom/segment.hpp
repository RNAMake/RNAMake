//
// Created by Joe Yesselman on 6/27/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_ALL_ATOM_SEGMENT_HPP_
#define RNAMAKE_SRC_STRUCTURE_ALL_ATOM_SEGMENT_HPP_

#include <fstream>

#include <structure/all_atom/basepair.h>
#include <structure/all_atom/residue.h>
#include <structure/base/segment.hpp>
#include <structure/secondary_structure/segment.hpp>
#include <structure/state/segment.hpp>

namespace structure::all_atom {
// typedefs to bring base into all_atom namespace
typedef structure::base::Chain<Residue> Chain;
typedef structure::base::Structure<Chain, Residue> Structure;
typedef structure::base::Pose<Basepair, Structure, Chain, Residue> Pose;
typedef structure::base::Segment<Basepair, Structure, Chain, Residue> Segment;

Segment get_segment_from_str(const String &);

secondary_structure::Segment get_secondary_structure(const Segment &);

state::Segment get_state(const Segment &);

void write_segment_to_pdb(const String &, const Segment &);

void align_segment(const Segment &, Segment &, Index);

} // namespace structure::all_atom

#endif // RNAMAKE_SRC_STRUCTURE_ALL_ATOM_SEGMENT_HPP_
