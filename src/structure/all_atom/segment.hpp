//
// Created by Joe Yesselman on 6/27/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_ALL_ATOM_SEGMENT_HPP_
#define RNAMAKE_SRC_STRUCTURE_ALL_ATOM_SEGMENT_HPP_

#include <structure/all_atom/basepair.h>
#include <structure/all_atom/residue.h>

namespace structure::all_atom {
// typedefs to bring base into all_atom namespace
typedef structure::Chain<Residue> Chain;
typedef structure::Structure<Chain, Residue> Structure;
typedef structure::Pose<Basepair, Structure, Chain, Residue> Pose;
typedef structure::Segment<Basepair, Structure, Chain, Residue> Segment;

Segment get_segment_from_str(const String &str) {
  Strings spl = base::string::split(str, "&");
  String name = spl[1];
  int aligned_end = std::stoi(spl[3]);
  auto mtype = static_cast<util::MotifType>(std::stoi(spl[4]));
  String structure_str = spl[5];
  Strings chain_strs = base::string::split(structure_str, ":");
  Residues res;
  structure::Cutpoints cutpoints;
  int i = 0;
  for (auto const &chain_str : chain_strs) {
    Strings res_strs = base::string::split(chain_str, ";");
    for (auto const &res_str : res_strs) {
      res.emplace_back(get_residue_from_str(res_str));
      i += 1;
    }
    cutpoints.push_back(i);
  }
  Structure s(res, cutpoints);
  Strings basepair_strs = base::string::split(spl[6], "@");
  Basepairs bps;
  for (auto const &bp_str : basepair_strs) {
    Strings bp_spl = base::string::split(bp_str, ",");
    Strings res_spl = base::string::split(bp_spl[0], "-");
    // TODO this does not account for chain ids now be longer than one char
    // We need to think of another strategy.
    String res1_id = res_spl[0].substr(0, 1);
    String res2_id = res_spl[1].substr(0, 1);
    int res1_num = std::stoi(res_spl[0].substr(1));
    int res2_num = std::stoi(res_spl[1].substr(1));
    const Residue &res1 = s.get_residue(res1_num, res1_id, ' ');
    const Residue &res2 = s.get_residue(res2_num, res2_id, ' ');
    Strings bp_state_strs = base::string::split(bp_spl[1], ";");
    math::Vector3 center = math::vector_from_str(bp_state_strs[0]);
    math::Matrix3x3 ref_frame = math::matrix_from_str(bp_state_strs[1]);
    math::Vector3s sugars;
    math::vectors_from_str(bp_state_strs[2], sugars);
    String bp_name = structure::generate_bp_name<Residue>(res1, res2);
    Basepair bp = {res1.get_uuid(),
                   res2.get_uuid(),
                   util::generate_uuid(),
                   structure::BasepairType::WC,
                   util::x3dna::X3dnaBPType::cWUW,
                   bp_name,
                   center,
                   sugars,
                   ref_frame};
    bps.emplace_back(bp);
  }
  Strings end_index_strs = base::string::split(spl[7], " ");
  Indexes end_indexes;
  for (auto const &strs : end_index_strs) {
    end_indexes.emplace_back(std::stoi(strs));
  }
  Strings end_ids = base::string::split(spl[8], ";");
  Strings ss_motif_str = base::string::split(spl[9], "!");
  String dot_bracket;
  for (auto const &chain_str : base::string::split(ss_motif_str[3], "|")) {
    for (auto const &res_str : base::string::split(chain_str, ";")) {
      Strings res_spl = base::string::split(res_str, ",");
      dot_bracket += res_spl[1];
    }
    dot_bracket += "&";
  }
  dot_bracket.pop_back();
  Structure proteins, small_molecules;
  // Pose p = {s, proteins, small_molecules, bps, end_indexes, end_ids, name};
  return {
      s,       proteins, small_molecules, bps,         end_indexes,
      end_ids, name,     mtype,           aligned_end, util::generate_uuid()
  };
}

} // namespace structure::all_atom

#endif // RNAMAKE_SRC_STRUCTURE_ALL_ATOM_SEGMENT_HPP_
