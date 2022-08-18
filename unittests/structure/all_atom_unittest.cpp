//
// Created by Joe Yesselman on 6/25/22.
//
#include "../common.hpp"

#include <base/paths.hpp>
#include <math/rotation.hpp>
#include <structure/all_atom/segment.hpp>

using namespace structure::all_atom;

/*
// TODO list of needed unittests
 - align(const Segment, Segment, index)

 // basepair fxns
 // TODO should we write a fxn for get_basepair_from_string() ?
 - new_uuids(r1 uuid, r2 uuid)
 - generate_bp_type(res1, res2, x3dna_bp_type)

 */

TEST_CASE("test all atom ") {
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
    structure::base::Cutpoints cutpoints;
    Structure s(res, cutpoints);
  }

  SUBCASE("test segment") {
    String path = base::path::resources_path() + "motifs/ref.motif";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    Segment seg = get_segment_from_str(lines[0]);
    structure::secondary_structure::Segment ss_seg =
        get_secondary_structure(seg);
    structure::state::Segment s_seg = get_state(seg);
  }
  SUBCASE("test alignment") {
    String path = base::path::resources_path() + "motifs/base.motif";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    Segment seg = get_segment_from_str(lines[0]);
    Segment seg2 = seg;
    math::Matrix3x3 rot;
    math::rotation_between_frames(seg.get_end_ref_frame(1),
                                  seg2.get_end_ref_frame(0), rot);
    seg2.rotate(rot);
    seg2.move(seg.get_end_center(1) - seg2.get_end_center(0));
    Real diff = math::difference_between_frames(seg.get_end_ref_frame(1),
                                                seg2.get_end_ref_frame(0));
    write_segment_to_pdb("test.pdb", seg);
    write_segment_to_pdb("test2.pdb", seg2);
  }
  SUBCASE("test alignment chain") {
    String path = base::path::resources_path() + "motifs/base.motif";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    std::vector<Segment> segs;
    segs.emplace_back(get_segment_from_str(lines[0]));
    for (int i = 0; i < 10; i++) {
      Segment new_seg = segs.back();
      align_segment(segs.back(), new_seg, 1);
      segs.emplace_back(new_seg);
    }
    for (int i = 0; i < segs.size(); i++) {
      write_segment_to_pdb("test." + std::to_string(i) + ".pdb", segs[i]);
    }
  }

  SUBCASE("test basepair") {
    String path = base::path::unittest_resource_path() +
                  "residue/test_str_to_residue.dat";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    // Residues res;

    SUBCASE("test basepair move") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      math::Vector3 vector_1 = {4, -1, 2};

      // TODO test move fxn
      // the two residue above should make a basepair
      // then the basepair should be moved by the specified vector
      // finally we should check the position of the atoms in the basepair
    }
    SUBCASE("test basepair rotate") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      math::Matrix3x3 matrix_1 = {-1, 0, 0, 0, 1, 0, 0, 0, 1};

      // TODO test rotate fxn
      // the two residue above should make a basepair
      // then the basepair should be rotated by the specified matrix
      // finally we should check the position of the atoms in the basepair
    }
    SUBCASE("test get uuid") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue
      util::Uuid uuid_1 = residue_2.get_uuid();

      // TODO test get_uuid fxn
      // that should be a basepair
      // find out what the UUID is and test it
    }
    SUBCASE("test get name") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue
      auto basepair_name = residue_2.get_name();

      // TODO test get_name fxn
      // should be a basepair
      // find out what the name is and test it
    }
    SUBCASE("test get center") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue
      auto center = residue_2.get_center();

      // TODO test get_center fxn
      // CHECK(center.get_x() == doctest::Approx(x));
      // CHECK(center.get_y() == doctest::Approx(y));
      // CHECK(center.get_z() == doctest::Approx(z));
    }
    SUBCASE("test get partner") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue

      // TODO test get_partner fxn
      // auto partner_uuid = basepair.get_partner(residue_1.get_uuid());
      // CHECK(partner_uuid == residue_2.get_uuid());
    }
    SUBCASE("test get basepair type") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue

      // TODO test get_bp_type fxn
      // should be a basepair
      // find out what the basepair type is and test it
    }
    SUBCASE("test swap residue positions") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue

      // TODO test swap_residue_positions
    }
    SUBCASE("test invert reference frame") {
      Residue residue_1 = get_residue_from_str(lines[0]);
      Residue residue_2 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue

      // TODO test invert_reference_frame
    }
    SUBCASE("test basepair constructors") {
      // TODO this should be a basepair

      // TODO test basepair constructors
    }
    SUBCASE("test is_equal") {
      // TODO test is_equal
      //Basepair basepair_1 = get_residue_from_str(lines[0]);
    }
    SUBCASE("test get res1 uuid") {
      // TODO test get_res_1_uuid
    }
    SUBCASE("test get res2 uuid") {
      // TODO test get_res_2_uuid
    }
    SUBCASE("test ref frame") {
      // TODO test get_ref_frame
    }
    SUBCASE("get res1_c1_prime_coord") {
      // TODO test get res1 c1 prime coord
    }
    SUBCASE("get res2_c1_prime_coord") {
      // TODO test get res2 c1 prime coord
    }
    SUBCASE("get c1_prime_coords") {
      // TODO test c1 prime coords
    }
    SUBCASE("") {

    }


  }
}