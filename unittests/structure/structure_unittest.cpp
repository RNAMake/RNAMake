//
// Created by Joe Yesselman on 6/25/22.
//
#include "../common.hpp"
#include <../structure/base/structure.hpp>
#include <../structure/secondary_structure/segment.hpp>

/*
// TODO list of tests needed here (from structure.hpp):
 - move(vector)
 - rotate(matrix)
*/

TEST_CASE("brief test of primitive functionality ") {
  SUBCASE("test constructor") {
    // TODO test constructor
  }
  SUBCASE("test begin") {
    // TODO test begin
  }
  SUBCASE("test get_residue with uuid") {
    // TODO test get_residue(uuid)
  }
  SUBCASE("test get_residue with index") {
    // TODO test get_residue(Index)
  }
  SUBCASE("test get_res_index with residue") {
    // TODO test get_res_index(residue)
  }
  SUBCASE("test get chains") {
    // TODO get_chains()
  }
  SUBCASE("test get cutpoints") {
    // TODO get_cutpoints()
  }
  SUBCASE("test get_num_residues()") {
    // TODO get_num_residues()
  }
  SUBCASE("test get_num_chains()") {
    // TODO get_num_chains()
  }
  SUBCASE("test get_sequence()") {
    // TODO get_sequence()
  }
  SUBCASE("test is residue start of chain") {
    // TODO test is_residue_start_of_chain()
  }
  SUBCASE("test is residue end of chain") {
    // TODO test is_residue end_of_chain()
  }
  SUBCASE("test move structure by vector") {
    // TODO test move(vector)
  }
  SUBCASE("test rotate structure by matrix") {
    // TODO test rotate(matrix)
  }

}