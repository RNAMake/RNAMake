//
// Created by Joseph Yesselman on 2019-04-11.
//

#ifndef RNAMAKE_NEW_SEQUENCE_TOOLS_H
#define RNAMAKE_NEW_SEQUENCE_TOOLS_H

#include <base/types.h>
#include <secondary_structure/pose.h>

namespace secondary_structure {

void
get_res_types_from_sequence(
        String const &,
        ResTypes & /* return array */);

int
find_res_types_in_pose(
        PoseOP,
        ResTypes const &);

int
find_gc_helix_stretches(
        PoseOP,
        int);

}

#endif //RNAMAKE_NEW_SEQUENCE_TOOLS_H
