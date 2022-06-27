//
// Created by Joe Yesselman on 6/25/22.
//

#include "../common.hpp"

#include <base/paths.hpp>
#include <structure/all_atom/basepair.h>
#include <structure/all_atom/residue.h>
#include <structure/base.hpp>

using namespace structure::all_atom;
Residue get_residue_from_str(String const &s) {
  Strings spl = base::string::split(s, ",");
  char name = spl[1][0];
  int num = std::stoi(spl[2]);
  String chain_id = spl[3];
  char i_code = ' ';
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
typedef structure::Pose<Basepair, Structure, Chain, Residue> Pose;
typedef structure::Segment<Basepair, Structure, Chain, Residue> Segment;

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
  SUBCASE("test basepair") {

  }
  SUBCASE("test segment") {
    String path = base::path::resources_path() + "motifs/ref.motif";
    auto lines = Strings();
    base::path::get_lines_from_file(path, lines);
    Strings spl = base::string::split(lines[0], "&");
    String m_path = spl[0];
    String name = spl[1];
    float score = std::stof(spl[2]);
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
      Basepair bp = {
          res1.get_uuid(),
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
    for(auto const & str : end_index_strs) {
      end_indexes.emplace_back(std::stoi(str));
    }
    Strings end_ids = base::string::split(spl[8], ";");
    Strings ss_motif_str = base::string::split(spl[9], "!");
    String dot_bracket;
    for(auto const & chain_str : base::string::split(ss_motif_str[3], "|")) {
      for(auto const & res_str : base::string::split(chain_str, ";")) {
        Strings res_spl = base::string::split(res_str, ",");
        dot_bracket += res_spl[1];
      }
      dot_bracket += "&";
    }
    dot_bracket.pop_back();
    Structure proteins, small_molecules;
    //Pose p = {s, proteins, small_molecules, bps, end_indexes, end_ids, name};
    Segment seg = {s, proteins, small_molecules, bps, end_indexes, end_ids,
        name, mtype, aligned_end, util::generate_uuid()};

  }
}