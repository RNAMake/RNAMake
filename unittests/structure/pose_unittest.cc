//
// Created by Joseph Yesselman on 12/3/17.
//

#include <iostream>
#include "../common.hpp"

#include <math/numerical.h>
#include <util/random_number_generator.h>
#include <structure/residue_type_set.h>
#include <structure/pdb_parser.h>
#include <structure/basepair.h>
#include <structure/pose.h>

TEST_CASE( "Test all atom pose") {
    auto rts = structure::ResidueTypeSet();
    auto parser = structure::PDBParser(rts);

    auto path = base::unittest_resources_path() + "/structure/p4p6.pdb";
    auto p = structure::get_pose_from_pdb(path, parser);


}
