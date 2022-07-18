//
// Created by Joe Yesselman on 7/6/22.
//

#include <structure/all_atom/segment.hpp>

namespace structure::all_atom {

Segment get_segment_from_str(const String &str) {
  Strings spl = ::base::string::split(str, "&");
  String name = spl[1];
  int aligned_end = std::stoi(spl[3]);
  auto mtype = static_cast<util::MotifType>(std::stoi(spl[4]));
  String structure_str = spl[5];
  Strings chain_strs = ::base::string::split(structure_str, ":");
  Residues res;
  structure::base::Cutpoints cutpoints;
  int i = 0;
  for (auto const &chain_str : chain_strs) {
    Strings res_strs = ::base::string::split(chain_str, ";");
    for (auto const &res_str : res_strs) {
      res.emplace_back(get_residue_from_str(res_str));
      i += 1;
    }
    cutpoints.push_back(res.size());
  }
  cutpoints.pop_back();
  Structure s(res, cutpoints);
  Strings basepair_strs = ::base::string::split(spl[6], "@");
  Basepairs bps;
  for (auto const &bp_str : basepair_strs) {
    Strings bp_spl = ::base::string::split(bp_str, ",");
    Strings res_spl = ::base::string::split(bp_spl[0], "-");
    // TODO this does not account for chain ids now be longer than one char
    // We need to think of another strategy.
    String res1_id = res_spl[0].substr(0, 1);
    String res2_id = res_spl[1].substr(0, 1);
    int res1_num = std::stoi(res_spl[0].substr(1));
    int res2_num = std::stoi(res_spl[1].substr(1));
    const Residue &res1 = s.get_residue(res1_num, res1_id, ' ');
    const Residue &res2 = s.get_residue(res2_num, res2_id, ' ');
    Strings bp_state_strs = ::base::string::split(bp_spl[1], ";");
    math::Vector3 center = math::vector_from_str(bp_state_strs[0]);
    int count = 0;
    // center computed wrong sometimes ... 
    math::Vector3 calc_center = {0, 0, 0};
    for(const auto & a : res1) {
      calc_center += a.get_coords();
      count += 1;
    }
    for(const auto & a : res2) {
      calc_center += a.get_coords();
      count += 1;
    }
    calc_center /= count;
    math::Matrix3x3 ref_frame = math::matrix_from_str(bp_state_strs[1]);
    math::Vector3s sugars;
    math::vectors_from_str(bp_state_strs[2], sugars);
    String bp_name = structure::base::generate_bp_name<Residue>(res1, res2);
    Basepair bp = {res1.get_uuid(),
                   res2.get_uuid(),
                   util::generate_uuid(),
                   structure::base::BasepairType::WC,
                   util::x3dna::X3dnaBPType::cWUW,
                   bp_name,
                   calc_center,
                   sugars,
                   ref_frame};
    bps.emplace_back(bp);
  }
  Strings end_index_strs = ::base::string::split(spl[7], " ");
  Indexes end_indexes;
  for (auto const &strs : end_index_strs) {
    end_indexes.emplace_back(std::stoi(strs));
  }
  Strings end_ids = ::base::string::split(spl[8], ";");
  Strings ss_motif_str = ::base::string::split(spl[9], "!");
  String dot_bracket;
  for (auto const &chain_str : ::base::string::split(ss_motif_str[3], "|")) {
    for (auto const &res_str : ::base::string::split(chain_str, ";")) {
      Strings res_spl = ::base::string::split(res_str, ",");
      dot_bracket += res_spl[1];
    }
    dot_bracket += "&";
  }
  dot_bracket.pop_back();
  Structure proteins, small_molecules;
  // Pose p = {s, proteins, small_molecules, bps, end_indexes, end_ids, name};
  return {s,
          proteins,
          small_molecules,
          bps,
          end_indexes,
          end_ids,
          name,
          dot_bracket,
          mtype,
          aligned_end,
          util::generate_uuid()};
}

secondary_structure::Segment get_secondary_structure(const Segment &seg) {
  String dot_bracket = seg.get_dot_bracket();
  secondary_structure::Residues res;
  int i = 0;
  for (auto const &r : seg) {
    String chain_id = r.get_chain_id();
    res.push_back({r.get_name(), dot_bracket[i], r.get_num(), chain_id,
                   r.get_i_code(), r.get_uuid(), r.get_rtype()});
    i += 1;
    if (seg.is_residue_end_of_chain(r)) {
      i += 1;
    }
  }
  structure::base::Cutpoints cuts = seg.get_cutpoints();
  secondary_structure::Structure s(res, cuts);
  secondary_structure::Basepairs bps;
  std::for_each(seg.bps_begin(), seg.bps_end(), [&bps](const Basepair &bp) {
    bps.push_back({bp.get_res1_uuid(), bp.get_res2_uuid(), bp.get_uuid(),
                   bp.get_bp_type()});
  });
  secondary_structure::Structure proteins, small_molecules;
  Indexes end_indexes = seg.get_end_indexes();
  Strings end_ids = seg.get_end_ids();
  util::MotifType mtype = seg.get_segment_type();
  Index aligned_end = seg.get_aligned_end_index();
  String name = seg.get_name();
  return {
      s,    proteins,    small_molecules, bps,         end_indexes,   end_ids,
      name, dot_bracket, mtype,           aligned_end, seg.get_uuid()};
}

state::Segment get_state(const Segment &seg) {
  state::Residues res;
  util::Beads beads;
  for (auto const &r : seg) {
    beads = r.get_beads();
    state::Residue s_r(beads);
    res.emplace_back(s_r);
  }
  structure::base::Cutpoints cuts = seg.get_cutpoints();
  state::Structure s(res, cuts);
  state::Basepairs bps;
  std::for_each(seg.bps_begin(), seg.bps_end(), [&bps](const Basepair &bp) {
    String name = bp.get_name();
    math::Vector3 center = bp.get_center();
    math::Vector3s c1_prime_coords = bp.get_c1_prime_coords();
    math::Matrix3x3 ref_frame = bp.get_ref_frame();
    bps.push_back({name, center, c1_prime_coords, ref_frame});
  });
  state::Structure proteins, small_molecules;
  String dot_bracket = seg.get_dot_bracket();
  Indexes end_indexes = seg.get_end_indexes();
  Strings end_ids = seg.get_end_ids();
  util::MotifType mtype = seg.get_segment_type();
  Index aligned_end = seg.get_aligned_end_index();
  String name = seg.get_name();
  return {
      s,    proteins,    small_molecules, bps,         end_indexes,   end_ids,
      name, dot_bracket, mtype,           aligned_end, seg.get_uuid()};
}

void write_segment_to_pdb(const String &fname, const Segment &seg) {
  std::ofstream out;
  out.open(fname);
  int acount = 1;
  int rnum = 1;
  char chain_id = 'A';
  for (auto const &r : seg) {
    for (auto const &a : r) {
      char buffer[200];
      math::Vector3 c = a.get_coords();
      std::sprintf(
          buffer,
          "%-6s%5d %-4s%1s%-4c%1c%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f     "
          " %4s%2s\n",
          "ATOM", acount, a.get_name().c_str(), "", r.get_name(), chain_id,
          rnum, "", c.get_x(), c.get_y(), c.get_z(), 1.00, 0.00, "", "");
      out << buffer;
      acount += 1;
    }
    rnum += 1;
  }
  out.close();
}

void align_segment(const Segment &ref, Segment &seg, Index end_index) {
  std::cout << seg.get_end_center(0) << std::endl;
  math::Matrix3x3 rot = math::rotation_between_frames(
      ref.get_end_ref_frame(end_index), seg.get_end_ref_frame(0));
  seg.rotate(rot);
  std::cout << seg.get_end_center(0) << std::endl;
  seg.move(ref.get_end_center(end_index) - seg.get_end_center(0));

}

}