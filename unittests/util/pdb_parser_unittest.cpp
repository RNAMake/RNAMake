//
// Created by Joe Yesselman on 6/23/22.
//

#include "../common.hpp"

#include <base/paths.hpp>
#include <util/io/pdb_parser.hpp>

TEST_CASE("test pdb parser") {
  using namespace util::io;
  SUBCASE("test trival") {
    String path =
        base::path::unittest_resource_path() + "/util/HELIX.IDEAL.pdb";
    std::cout << path << std::endl;
    PDBParser pdb_parser;
    pdb_parser.parse(path);
    CHECK(pdb_parser.done() == false);
    pdb_parser.next();
    CHECK(pdb_parser.get_line_type() == PDBLineType::ATOM);
    auto const & atom_data = pdb_parser.get_line_atom_data();
    CHECK(atom_data.atom_id == 1);
    CHECK(atom_data.atom_name == "P");
    CHECK(atom_data.res_name == "G");
    CHECK(atom_data.res_num == 4);
    CHECK(atom_data.x == doctest::Approx(7.719));
    CHECK(atom_data.y == doctest::Approx(4.035));
    CHECK(atom_data.z == doctest::Approx(12.18));
    int count = 0;
    while(!pdb_parser.done()) {
      pdb_parser.next();
      count += 1;
    }
    CHECK(count == 88);
  }
  SUBCASE("test improper pdb format") {
    
  }
}