//
// Created by Joe Yesselman on 6/17/22.
//

#include "../common.hpp"
#include <util/find_pair.h>
#include <util/x3dna/x3dna.h>

TEST_CASE("find pair unittests ") {
  SUBCASE("test trivial") {
    auto os_name = "osx";
    auto x3dna_path = base::path::resources_path() + "/x3dna/" + os_name + "/";
    String env = "X3DNA=" + x3dna_path;
    auto s = strdup(env.c_str());
    putenv(s);
    util::PairFinder pf("/Users/jyesselm/projects/RNAMake/unittests/unittest_resources/util/HELIX.IDEAL.pdb");
    util::x3dna::X3dna::X3Basepairs pairs;
    pf.find_pair(pairs);
    CHECK(pairs.size() == doctest::Approx(2));
  }
}
