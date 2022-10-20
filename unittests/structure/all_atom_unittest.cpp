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
    Residue residue_1 = get_residue_from_str(lines[0]);
    Residue residue_2 = get_residue_from_str(lines[1]);
    math::Vector3 res_1_center = residue_1.get_center();
    math::Vector3 res_2_center = residue_2.get_center();
    math::Vector3 basepair_1_center;
    basepair_1_center.set_x((res_1_center.get_x() + res_2_center.get_x()) / 2);
    basepair_1_center.set_y((res_1_center.get_y() + res_2_center.get_y()) / 2);
    basepair_1_center.set_z((res_1_center.get_z() + res_2_center.get_z()) / 2);
    math::Matrix3x3 basepair_ref_frame;
    structure::base::BasepairType basepair_type =
        structure::base::BasepairType::WC;
    auto basepair_1_name = "test basepair";
    auto basepair_x3dna_type = util::x3dna::get_x3dna_by_type("cM-");
    // TODO these are not the actual C1' coords, just filler for debugging
    auto c1_prime_coords = {res_1_center, res_2_center};

    Basepair basepair_1 = Basepair(
        residue_1.get_uuid(), residue_2.get_uuid(), util::generate_uuid(),
        basepair_type, basepair_x3dna_type, basepair_1_name, basepair_1_center,
        c1_prime_coords, basepair_ref_frame);

    SUBCASE("test basepair move") {
      math::Vector3 vector_1 = {4, -1, 2};
      basepair_1.move(vector_1);
      basepair_1_center = basepair_1.get_center();

      CHECK(basepair_1_center.get_x() == doctest::Approx(4));
      CHECK(basepair_1_center.get_y() == doctest::Approx(-1));
      CHECK(basepair_1_center.get_z() == doctest::Approx(2));

      math::Vector3 vector_2 = {-4, 1, -2};
      basepair_1.move(vector_1);
      basepair_1_center = basepair_1.get_center();

      CHECK(basepair_1_center.get_x() == doctest::Approx(0));
      CHECK(basepair_1_center.get_y() == doctest::Approx(0));
      CHECK(basepair_1_center.get_z() == doctest::Approx(0));
    }
    SUBCASE("test basepair rotate") {
      SUBCASE("test basepair rotate x") {
        math::Matrix3x3 matrix_1 = {-1, 0, 0, 0, 1, 0, 0, 0, 1};
        basepair_1.rotate(matrix_1);
        math::Matrix3x3 ref_frame_matrix = basepair_1.get_ref_frame();

        CHECK(ref_frame_matrix.get_xx() == doctest::Approx(-1));
        CHECK(ref_frame_matrix.get_xy() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_xz() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_yx() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_yy() == doctest::Approx(1));
        CHECK(ref_frame_matrix.get_yz() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zx() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zy() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zz() == doctest::Approx(1));
      }
      SUBCASE("test basepair rotate y") {
        math::Matrix3x3 matrix_1 = {1, 0, 0, 0, -1, 0, 0, 0, 1};
        basepair_1.rotate(matrix_1);
        math::Matrix3x3 ref_frame_matrix = basepair_1.get_ref_frame();

        CHECK(ref_frame_matrix.get_xx() == doctest::Approx(1));
        CHECK(ref_frame_matrix.get_xy() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_xz() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_yx() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_yy() == doctest::Approx(-1));
        CHECK(ref_frame_matrix.get_yz() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zx() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zy() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zz() == doctest::Approx(1));
      }
      SUBCASE("test basepair rotate z") {
        math::Matrix3x3 matrix_1 = {1, 0, 0, 0, 1, 0, 0, 0, -1};
        basepair_1.rotate(matrix_1);
        math::Matrix3x3 ref_frame_matrix = basepair_1.get_ref_frame();

        CHECK(ref_frame_matrix.get_xx() == doctest::Approx(1));
        CHECK(ref_frame_matrix.get_xy() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_xz() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_yx() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_yy() == doctest::Approx(1));
        CHECK(ref_frame_matrix.get_yz() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zx() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zy() == doctest::Approx(0));
        CHECK(ref_frame_matrix.get_zz() == doctest::Approx(-1));
      }
    }
    SUBCASE("test get uuid") {
      Residue residue_3 = get_residue_from_str(lines[0]);
      Residue residue_4 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue
      util::Uuid uuid_1 = residue_2.get_uuid();

      // TODO test get_uuid fxn
      // that should be a basepair
      // find out what the UUID is and test it
    }
    SUBCASE("test get name") {
      Residue residue_3 = get_residue_from_str(lines[0]);
      Residue residue_4 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue
      auto basepair_name = residue_2.get_name();

      // TODO test get_name fxn
      // should be a basepair
      // find out what the name is and test it
    }
    SUBCASE("test get center") {
      Residue residue_3 = get_residue_from_str(lines[0]);
      Residue residue_4 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue
      auto center = residue_2.get_center();

      // TODO test get_center fxn
      // CHECK(center.get_x() == doctest::Approx(x));
      // CHECK(center.get_y() == doctest::Approx(y));
      // CHECK(center.get_z() == doctest::Approx(z));
    }
    SUBCASE("test get partner") {
      Residue residue_3 = get_residue_from_str(lines[0]);
      Residue residue_4 = get_residue_from_str(lines[1]);
      // TODO this should be a basepair not a residue

      // TODO test get_partner fxn
      // auto partner_uuid = basepair.get_partner(residue_1.get_uuid());
      // CHECK(partner_uuid == residue_2.get_uuid());
    }
    SUBCASE("test get basepair type") {
      // TODO test get_bp_type fxn
      // should be a basepair
      // find out what the basepair type is and test it
    }
    SUBCASE("test swap residue positions") {
      // TODO test swap_residue_positions
    }
    SUBCASE("test invert reference frame") {
      math::Vector3 vector_1 = {4, -2, 1};
      basepair_1.move(vector_1);
      basepair_1.invert_reference_frame();
      math::Matrix3x3 basepair_1_ref_frame = basepair_1.get_ref_frame();

      CHECK(basepair_1_ref_frame.get_xx() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_xy() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_xz() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_yx() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_yy() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_yz() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_zx() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_zy() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_zz() == doctest::Approx(0));
      // TODO fix this so the reference frame is actually inverted
    }
    SUBCASE("test basepair constructors") {
      // TODO test basepair constructors
    }
    SUBCASE("test is_equal") {
      SUBCASE("test is_equal = true") {}
      SUBCASE("test is_equal = false") {}
    }
    SUBCASE("test get res1 uuid") {
      auto uuid_1 = basepair_1.get_res1_uuid();
      CHECK(basepair_1.get_res1_uuid() == uuid_1);
    }
    SUBCASE("test get res2 uuid") {
      auto uuid_2 = basepair_1.get_res2_uuid();
      CHECK(basepair_1.get_res2_uuid() == uuid_2);
    }
    SUBCASE("test get ref frame") {
      math::Matrix3x3 basepair_1_ref_frame = basepair_1.get_ref_frame();

      CHECK(basepair_1_ref_frame.get_xx() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_xy() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_xz() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_yx() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_yy() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_yz() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_zx() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_zy() == doctest::Approx(0));
      CHECK(basepair_1_ref_frame.get_zz() == doctest::Approx(0));
      // TODO double check to make sure the right ref_frame coords are in the
      // doctest expression
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
    SUBCASE("") {}
  }
}