//
// Created by Joe Yesselman on 6/25/22.
//

#include "../common.hpp"

#include <base/paths.hpp>
#include <structure/base.hpp>
#include <structure/all_atom/residue.h>

using namespace structure::all_atom;
Residue get_residue_from_str(String const &s) {
  Strings spl = base::string::split(s, ",");
  char name = spl[1][0];
  int num = std::stoi(spl[2]);
  String chain_id = spl[3];
  char i_code = spl[4][0];
  util::Uuid uuid = util::generate_uuid();
  Atoms atoms;
  int i = 5;
  while (i < spl.size()) {
    if (spl[i].size() <= 1) {
      i++;
      continue;
    }
    atoms.push_back(Atom(spl[i]));
    i++;
  }
  return {name, num, chain_id, i_code, atoms, uuid};
}

typedef structure::Chain<Residue> Chain;
typedef structure::Structure<Chain, Residue> Structure;

TEST_CASE("test all atom") {
  SUBCASE("test generate residue") {
    String path = base::path::unittest_resource_path() +
                  "residue/test_str_to_residue.dat";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    Residue r = get_residue_from_str(lines[0]);
  }
  SUBCASE("test chain") {
    String path = base::path::unittest_resource_path() +
                  "residue/test_str_to_residue.dat";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    Residues res;
    res.emplace_back(get_residue_from_str(lines[0]));
    res.emplace_back(get_residue_from_str(lines[1]));
    Chain c(res);
    CHECK(c.get_first().get_num_atoms() == res[0].get_num_atoms());
    CHECK(c.get_last().get_num_atoms() == res[1].get_num_atoms());
  }
  SUBCASE("test structure") {
    String path = base::path::unittest_resource_path() +
                  "residue/test_str_to_residue.dat";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    Residues res;
    res.emplace_back(get_residue_from_str(lines[0]));
    res.emplace_back(get_residue_from_str(lines[1]));
    structure::Cutpoints cutpoints;
    Structure s(res, cutpoints);
  }
}