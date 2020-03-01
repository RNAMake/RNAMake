//
// Created by Joseph Yesselman on 12/7/17.
//

#include <iostream>
#include "../common.hpp"

#include <base/settings.h>
#include <base/file_io.h>
#include <util/x3dna.h>

TEST_CASE( "Test X3dna parser ", "[X3dnaParser]" ) {
    //init_unittest_safe_logging();

    auto x = util::X3dna();
    auto path = base::unittest_resource_dir() + "/util/p4p6.pdb";
    auto x3dna_bps = x.get_basepairs(path);
    REQUIRE(x3dna_bps.size() == 63);
    REQUIRE(!base::file_exists("ref_frames.dat"));
    REQUIRE(!base::file_exists("p4p6_dssr.out"));

}