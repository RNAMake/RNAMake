//
// Created by Joseph Yesselman on 12/7/17.
//

#include <iostream>
#include "../common.hpp"

#include <base/settings.h>
#include <base/file_io.h>
#include <util/x3dna.h>
#include <util/dssr.h>
#include <base/sys_interface.h>

TEST_CASE( "Test X3dna parser ", "[X3dnaParser]" ) {
    //init_unittest_safe_logging();

    auto x = util::X3dna();
    auto path = base::unittest_resource_dir() + "/util/p4p6.pdb";
    auto x3dna_bps = x.get_basepairs(path);
    auto nts = util::DssrNts{}; 
    auto bps = util::DssrPairs{}; 
    get_elements(path,nts,bps);
    REQUIRE(x3dna_bps.size() == 63);
    REQUIRE(!base::file_exists("ref_frames.dat"));
    REQUIRE(!base::file_exists("p4p6_dssr.out"));

}

TEST_CASE( "Test X3dna parser JSON version", "[X3dnaParser]" ) {
    //init_unittest_safe_logging();
    auto x = util::X3dna();
    auto path = base::unittest_resource_dir() + "/util/p4p6.pdb";
    //auto x3dna_bps = x.get_basepairs_json("../../data/255D.cif");
    auto x3dna_bps = x.get_basepairs_json(path);
    REQUIRE(x3dna_bps.size() == 76);
    REQUIRE(!base::file_exists("ref_frames.dat"));
    REQUIRE(!base::file_exists("p4p6_dssr.out"));
    

}

TEST_CASE( "Comparing outputs for some small structures", "[X3dnaParser]" ) {
    auto x = util::X3dna();
    for(auto pdb : {"../../data/100D.pdb","../../data/157D.pdb","../../data/1DI2.pdb"}){
        auto x3dna_bps = x.get_basepairs(pdb);
        auto x3dna_bps_json = x.get_basepairs_json(pdb);
        util::compare_bps(x3dna_bps,x3dna_bps_json); 
    }
    REQUIRE(true);

}
