//
// Created by Joseph Yesselman on 12/7/17.
//

#include <iostream>
#include "../common.hpp"

#include <base/settings.h>
#include <base/file_io.h>
#include <util/x3dna.h>
#include <util/dssr.h>
#include "base/env_manager.h"
#include <base/sys_interface.h>

/* so what are the things we want to test? 
 *
 *
 */

TEST_CASE("Testing basepair code converters ", "[get_x3dna_by_type,get_str_from_x3dna_type]" ) {
  util::X3dna::set_envs();


  const auto codes = Strings{"cm-" ,"cM-M","tW+W","c.+M",".W+W","tW-M","tm-M","cW+M",".W-W","cM+.","c.-m","cM+W","tM+m","tM-W","cm-m","cM-W","cW-W","c.-M","cm+M","cm-M","....","cm-W","tM-m","c.-W","cM+m","cM-m","c...","tW+m","c.+m","tm+m","tW+.","tm+W","t...","cW-.","cW-M","t.-W","tM+M","t.-M","cM-.","cW-m","t.+m","tM-.","cm+W","cM+M","cm+.","cm-.","c.-.","cW+W","t.-.","t.+W","tm-m","cW+.","tm+.","t.+.","c.+.","t.-m","t.+M","tW-.","tm-W","tM-M","tM+.","c.+W","tm+M","tW-m","cW+m","tm-.","tW+M",".W+m","tM+W","..+m","tW-W","cm+m",".W-m",".M+m",".W+M",".M+M",".m+W",".W-M",".m+m","..-M",".M-m","..-m",".M+.",".m-m",".M-W",".W-."};


    for(const auto& bp_code : codes) {
        REQUIRE(bp_code == util::get_str_from_x3dna_type(util::get_x3dna_by_type(bp_code)));
    }
}
 TEST_CASE( "Test X3dna parser ", "[X3dnaParser]" ) {
    util::X3dna::set_envs();

  //init_unittest_safe_logging();

    auto x = util::X3dna();
    auto path = base::unittest_resource_dir() + "/util/p4p6.pdb";
    auto x3dna_bps = x.get_basepairs(path);
    auto nts = util::DssrNts{}; 
    auto bps = util::DssrPairs{}; 
    auto hps = util::DssrHairpins{}; 
    auto hels = util::DssrHelices{}; 
    auto stems = util::DssrStems{}; 
    auto iloops = util::DssrILoops{}; 
    get_elements(path,nts,bps,hps,hels,stems,iloops);
    REQUIRE(x3dna_bps.size() == 63);
    REQUIRE(!base::file_exists("ref_frames.dat"));
    REQUIRE(!base::file_exists("p4p6_dssr.out"));

}

TEST_CASE( "Test X3dna parser JSON version", "[X3dnaParser]" ) {
    //init_unittest_safe_logging();
    auto x = util::X3dna();
    auto path = base::unittest_resource_dir() + "/util/p4p6.pdb";
    auto x3dna_bps = x.get_basepairs_json(path);
    REQUIRE(x3dna_bps.size() == 76);
    REQUIRE(!base::file_exists("ref_frames.dat"));
    REQUIRE(!base::file_exists("p4p6_dssr.out"));
    

}


