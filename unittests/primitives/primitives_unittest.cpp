//
// Created by Joe Yesselman on 6/22/22.
//

#include "../common.hpp"

#include <primitives/basepair.h>
#include <primitives/residue.h>

TEST_CASE("brief test of primitive functionility") {
  SUBCASE("test residue") {
    primitives::Residue r1 = {'A', 1, 'A', ' ', util::generate_uuid()};
    CHECK(r1.get_name() == 'A');
    CHECK(r1.get_num() == 1);
    auto uuid = r1.get_uuid();
    CHECK(r1.get_uuid() == uuid);
    primitives::Residue r2 = {'A', 1, 'A', ' ', uuid};
    CHECK(r1 == r2);
  }
  SUBCASE("test basepairs") {
    SUBCASE("test default construction") {
      String name = "A1-A2";
      primitives::Basepair bp1 = {util::generate_uuid(), util::generate_uuid(),
                                  util::generate_uuid(),
                                  primitives::BasepairType::WC, name};
      CHECK(bp1.get_name() == "A1-A2");
    }
    SUBCASE("test name construction") {
      primitives::Residue r1 = {'A', 1, 'A', ' ', util::generate_uuid()};
      
    }
  }
}