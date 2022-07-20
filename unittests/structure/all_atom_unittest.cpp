//
// Created by Joe Yesselman on 6/25/22.
//
#include "../common.hpp"

#include <base/paths.hpp>
#include <math/rotation.hpp>
#include <structure/all_atom/segment.hpp>

using namespace structure::all_atom;

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
  SUBCASE("test basepair") {}
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
    // std::cout << seg.get_end(1).get_ref_frame().get_flip_orientation()
    //                  .get_str() << std::endl;
    math::rotation_between_frames(seg.get_end_ref_frame(1),
                                  seg2.get_end_ref_frame(0), rot);
    seg2.rotate(rot);
    seg2.move(seg.get_end_center(1) - seg2.get_end_center(0));
    Real diff = math::difference_between_frames(seg.get_end_ref_frame(1),
                                                seg2.get_end_ref_frame(0));
    //std::cout << diff*180.0 / PI << std::endl;
    write_segment_to_pdb("test.pdb", seg);
    write_segment_to_pdb("test2.pdb", seg2);
  }
  SUBCASE("test alignment chain") {
    String path = base::path::resources_path() + "motifs/base.motif";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    std::vector<Segment> segs;
    segs.emplace_back(get_segment_from_str(lines[0]));
    for(int i = 0; i < 10; i++) {
      Segment new_seg = segs.back();
      align_segment(segs.back(), new_seg, 1);
      segs.emplace_back(new_seg);
    }
    for(int i =0 ; i < segs.size(); i++) {
      write_segment_to_pdb("test."+std::to_string(i)+".pdb", segs[i]);
    }
  }
}