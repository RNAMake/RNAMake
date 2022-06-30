//
// Created by Joe Yesselman on 6/30/22.
//
#include "../common.hpp"

#include <base/paths.hpp>
#include <math/rotation.hpp>
#include <structure/all_atom/residue.h>

using namespace structure::all_atom;
TEST_CASE("test all atom") {
  /*SUBCASE("test string conversion consistency") {
    for(int i =0 ; i < 1; i++) {
      String path = base::path::unittest_resource_path() +
                    "residue/test_str_to_residue.dat";
      auto lines = Strings();
      base::path::get_lines_from_file(path, lines);
      Residue r = get_residue_from_str(lines[0]);
      math::Vector3 v = {rand(), rand(), rand()};
      r.move(v);
    }
  } */
  SUBCASE("test string conversion consistency from seg str") {
    for(int a = 0; a < 1000; a++) {
      String path = base::path::resources_path() + "motifs/base.motif";
      auto lines = Strings();
      base::path::get_lines_from_file(path, lines);
      Strings spl = base::string::split(lines[0], "&");
      std::cout << "ORG STRING:" << lines[0] << std::endl;
      Strings chain_strs = base::string::split(spl[5], ":");
      //CHECK(chain_strs[0].length() == 3250);
      //CHECK(chain_strs[1].length() == 2708);
      Residues res;
      for (auto const &chain_str : chain_strs) {
        Strings res_strs = base::string::split(chain_str, ";");
        for (auto const &res_str : res_strs) {
          res.push_back(get_residue_from_str(res_str));
        }
      }
      math::Vector3 v = {rand(), rand(), rand()};
      res[0].move(v);
    }
  }

}
