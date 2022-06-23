//
// Created by Joe Yesselman on 6/22/22.
//
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include <map>
#include <set>
#include <string>
#include <vector>

// RNAMake Headers
#include <math/matrix_3x3.hpp>
#include <util/x3dna/x3dna.h>

// contents:
// ana_fncs.cpp : line ~730
// app_fncs.cpp : line ~5980
// cmn_fncs.cpp : line ~8070
// analyze.cpp : line ~13740
// nrutil.cpp : line ~14500

using namespace std;

namespace util {
namespace x3dna {
X3dna::X3dna()
    : rebuild_files_(true), generated_dssr_(false),
      generated_ref_frames_(false), no_ref_frames_(false) {

  auto os_name = getprogname();
}

void X3dna::_delete_file(String const &file_name) const {
  try {
    std::remove(file_name.c_str());
  } catch (String const &e) {
  }
}

math::Vector3 X3dna::_convert_string_to_point(String const &str) const {
  auto doubles = std::vector<double>();
  auto spl = base::string::split(str, " ");
  for (auto const &s : spl) {
    if (s.length() > 1) {
      doubles.push_back(std::stod(s));
    }
    if (doubles.size() == 3) {
      break;
    }
  }
  return math::Vector3(doubles[0], doubles[1], doubles[2]);
}

void X3dna::_parse_ref_frame_file(String const &pdb_path,
                                  X3Basepairs &basepairs) const {
  // auto finder = PairFinder(pdb_path);
  // finder.find_pair(basepairs);
}

Strings X3dna::_split_over_white_space(String const &str) const {
  Strings spl = base::string::split(str, " ");
  Strings non_white_space;
  String temp;
  for (auto &s : spl) {
    temp = base::string::trim(s);
    if (temp.size() == 0) {
      continue;
    }
    non_white_space.push_back(temp);
  }
  return non_white_space;
}

/*
X3dna::X3Basepairs X3dna::get_basepairs_json(String const &pdb_path) const {

  auto dssr_json = base::execute_command_json(
      bin_path_ + "/x3dna-dssr -i=" + pdb_path + " --json --more 2>.error");
  // deleting the temp files that we don't want.
  util::json_cleanup();
  auto nt_it = dssr_json.find("nts");
  if (nt_it == dssr_json.end() || nt_it->is_null() || nt_it->empty()) {
    return X3Basepairs{};
  }
  // loop through the nucleotides and store them in a map
  auto res_map = std::map<String, X3Residue>{};
  for (auto &nt : *dssr_json.find("nts")) {
    auto residue_it = util::get_string(nt, "nt_id");
    auto ii = util::get_int(nt, "nt_resnum");
    auto chain = util::get_char(nt, "chain_name");
    res_map[residue_it] = X3Residue(ii, chain, ' ');
  }
  // now the basepairs can be constructed
  auto basepairs = X3Basepairs();
  auto pairings = dssr_json.find("pairs");
  if (pairings == dssr_json.end() || pairings->empty()) {
    return basepairs;
  }
  // loop through the pairings
  for (auto pp : *pairings) {
    auto name = util::get_string(pp, "name");
    auto ii = util::get_string(pp, "Saenger");
    auto nt_1 = util::get_string(pp, "nt1");
    auto nt_2 = util::get_string(pp, "nt2");
    auto type = util::get_string(pp, "DSSR");
    // get the reference frame information
    auto frame = math::Matrix3x3{};
    auto origin = math::Vector3{};
    auto frame_it = pp.find("frame");

    if (frame_it != pp.end() && !frame_it->is_null()) {
      origin = get_point(*frame_it, "origin");
      frame = get_matrix(*frame_it);
    }
    // define the basepair type and residues
    auto bp_type = get_x3dna_by_type(type);
    auto res1 = res_map.at(nt_1); // TODO maybe add an unkown thing?
    auto res2 = res_map.at(nt_2);
    // TODO add check that they are correct
    // add the new basepair i
    auto new_bp = X3Basepair{res1, res2, origin, frame, bp_type};
    // std::cout<<nt_1<<"\t"<<nt_2<<"\t"<<type<<std::endl;
    if (new_bp.valid()) {
      basepairs.push_back(std::move(new_bp));
    }
  }

  const auto num_basepairs = util::get_int(dssr_json, "num_pairs");
  if (num_basepairs != basepairs.size()) {
    throw X3dnaException("Basepair count mismatch. Expected " +
                         std::to_string(num_basepairs) + " but got " +
                         std::to_string(basepairs.size()));
  }
  return basepairs;
}
*/

X3dna::X3Basepairs X3dna::get_basepairs(String const &pdb_path) const {
  auto basepairs = X3Basepairs();
  _parse_ref_frame_file(pdb_path, basepairs);
  return basepairs;
}

X3dnaBPType get_x3dna_by_type(String const &name) {
  if (name == "cm-") {
    return X3dnaBPType::cmU;
  } else if (name == "cM-M") {
    return X3dnaBPType::cMUM;
  } else if (name == "tW+W") {
    return X3dnaBPType::tWPW;
  } else if (name == "c.+M") {
    return X3dnaBPType::cDPM;
  } else if (name == ".W+W") {
    return X3dnaBPType::DWPW;
  } else if (name == "tW-M") {
    return X3dnaBPType::tWUM;
  } else if (name == "tm-M") {
    return X3dnaBPType::tmUM;
  } else if (name == "cW+M") {
    return X3dnaBPType::cWPM;
  } else if (name == ".W-W") {
    return X3dnaBPType::DWUW;
  } else if (name == "cM+.") {
    return X3dnaBPType::cMPD;
  } else if (name == "c.-m") {
    return X3dnaBPType::cDUm;
  } else if (name == "cM+W") {
    return X3dnaBPType::cMPW;
  } else if (name == "tM+m") {
    return X3dnaBPType::tMPm;
  } else if (name == "tM-W") {
    return X3dnaBPType::tMUW;
  } else if (name == "cm-m") {
    return X3dnaBPType::cmUm;
  } else if (name == "cM-W") {
    return X3dnaBPType::cMUW;
  } else if (name == "cW-W") {
    return X3dnaBPType::cWUW;
  } else if (name == "c.-M") {
    return X3dnaBPType::cDUM;
  } else if (name == "cm+M") {
    return X3dnaBPType::cmPM;
  } else if (name == "cm-M") {
    return X3dnaBPType::cmUM;
  } else if (name == "....") {
    return X3dnaBPType::DDDD;
  } else if (name == "cm-W") {
    return X3dnaBPType::cmUW;
  } else if (name == "tM-m") {
    return X3dnaBPType::tMUm;
  } else if (name == "c.-W") {
    return X3dnaBPType::cDUW;
  } else if (name == "cM+m") {
    return X3dnaBPType::cMPm;
  } else if (name == "cM-m") {
    return X3dnaBPType::cMUm;
  } else if (name == "c...") {
    return X3dnaBPType::cDDD;
  } else if (name == "tW+m") {
    return X3dnaBPType::tWPm;
  } else if (name == "c.+m") {
    return X3dnaBPType::cDPm;
  } else if (name == "tm+m") {
    return X3dnaBPType::tmPm;
  } else if (name == "tW+.") {
    return X3dnaBPType::tWPD;
  } else if (name == "tm+W") {
    return X3dnaBPType::tmPW;
  } else if (name == "t...") {
    return X3dnaBPType::tDDD;
  } else if (name == "cW-.") {
    return X3dnaBPType::cWUD;
  } else if (name == "cW-M") {
    return X3dnaBPType::cWUM;
  } else if (name == "t.-W") {
    return X3dnaBPType::tDUW;
  } else if (name == "tM+M") {
    return X3dnaBPType::tMPM;
  } else if (name == "t.-M") {
    return X3dnaBPType::tDUM;
  } else if (name == "cM-.") {
    return X3dnaBPType::cMUD;
  } else if (name == "cW-m") {
    return X3dnaBPType::cWUm;
  } else if (name == "t.+m") {
    return X3dnaBPType::tDPm;
  } else if (name == "tM-.") {
    return X3dnaBPType::tMUD;
  } else if (name == "cm+W") {
    return X3dnaBPType::cmPW;
  } else if (name == "cM+M") {
    return X3dnaBPType::cMPM;
  } else if (name == "cm+.") {
    return X3dnaBPType::cmPD;
  } else if (name == "cm-.") {
    return X3dnaBPType::cmUD;
  } else if (name == "c.-.") {
    return X3dnaBPType::cDUD;
  } else if (name == "cW+W") {
    return X3dnaBPType::cWPW;
  } else if (name == "t.-.") {
    return X3dnaBPType::tDUD;
  } else if (name == "t.+W") {
    return X3dnaBPType::tDPW;
  } else if (name == "tm-m") {
    return X3dnaBPType::tmUm;
  } else if (name == "cW+.") {
    return X3dnaBPType::cWPD;
  } else if (name == "tm+.") {
    return X3dnaBPType::tmPD;
  } else if (name == "t.+.") {
    return X3dnaBPType::tDPD;
  } else if (name == "c.+.") {
    return X3dnaBPType::cDPD;
  } else if (name == "t.-m") {
    return X3dnaBPType::tDUm;
  } else if (name == "t.+M") {
    return X3dnaBPType::tDPM;
  }
  // added by CJ
  else if (name == "tW-.") {
    return X3dnaBPType::tWUD;
  } else if (name == "tm-W") {
    return X3dnaBPType::tmUW;
  } else if (name == "tM-M") {
    return X3dnaBPType::tMUM;
  } else if (name == "tM+.") {
    return X3dnaBPType::tMPD;
  } else if (name == "c.+W") {
    return X3dnaBPType::cDPW;
  } else if (name == "tm+M") {
    return X3dnaBPType::tmPM;
  } else if (name == "tW-m") {
    return X3dnaBPType::tWUm;
  } else if (name == "cW+m") {
    return X3dnaBPType::cWPm;
  } else if (name == "tm-.") {
    return X3dnaBPType::tmUD;
  } else if (name == "tW+M") {
    return X3dnaBPType::tWPM;
  } else if (name == ".W+m") {
    return X3dnaBPType::DWPm;
  } else if (name == "tM+W") {
    return X3dnaBPType::tMPW;
  } else if (name == "..+m") {
    return X3dnaBPType::DDPm;
  } else if (name == "tW-W") {
    return X3dnaBPType::tWUW;
  } else if (name == "cm+m") {
    return X3dnaBPType::cmPm;
  } else if (name == ".W-m") {
    return X3dnaBPType::DWUm;
  } else if (name == ".M+m") {
    return X3dnaBPType::DMPm;
  } else if (name == ".W+M") {
    return X3dnaBPType::DWPM;
  } else if (name == ".M+M") {
    return X3dnaBPType::DMPM;
  } else if (name == ".m+W") {
    return X3dnaBPType::DmPW;
  } else if (name == ".W-M") {
    return X3dnaBPType::DWUM;
  } else if (name == ".m+m") {
    return X3dnaBPType::DmPm;
  } else if (name == "..-M") {
    return X3dnaBPType::DDUM;
  } else if (name == ".M-m") {
    return X3dnaBPType::DMUm;
  } else if (name == "..-m") {
    return X3dnaBPType::DDUm;
  } else if (name == ".M+W") {
    return X3dnaBPType::DMPW;
  } else if (name == ".M+.") {
    return X3dnaBPType::DMPD;
  } else if (name == ".M-M") {
    return X3dnaBPType::DMUM;
  } else if (name == ".m-m") {
    return X3dnaBPType::DmUm;
  } else if (name == ".M-W") {
    return X3dnaBPType::DMUW;
  } else if (name == ".W-.") {
    return X3dnaBPType::DWUD;
  } else {
    throw X3dnaException("cannot get x3dna type with: " + name);
  }
}

String get_str_from_x3dna_type(X3dnaBPType type) {
  if (type == X3dnaBPType::cmU) {
    return "cm-";
  } else if (type == X3dnaBPType::cMUM) {
    return "cM-M";
  } else if (type == X3dnaBPType::tWPW) {
    return "tW+W";
  } else if (type == X3dnaBPType::cDPM) {
    return "c.+M";
  } else if (type == X3dnaBPType::DWPW) {
    return ".W+W";
  } else if (type == X3dnaBPType::tWUM) {
    return "tW-M";
  } else if (type == X3dnaBPType::tmUM) {
    return "tm-M";
  } else if (type == X3dnaBPType::cWPM) {
    return "cW+M";
  } else if (type == X3dnaBPType::DWUW) {
    return ".W-W";
  } else if (type == X3dnaBPType::cMPD) {
    return "cM+.";
  } else if (type == X3dnaBPType::cDUm) {
    return "c.-m";
  } else if (type == X3dnaBPType::cMPW) {
    return "cM+W";
  } else if (type == X3dnaBPType::tMPm) {
    return "tM+m";
  } else if (type == X3dnaBPType::tMUW) {
    return "tM-W";
  } else if (type == X3dnaBPType::cmUm) {
    return "cm-m";
  } else if (type == X3dnaBPType::cMUW) {
    return "cM-W";
  } else if (type == X3dnaBPType::cWUW) {
    return "cW-W";
  } else if (type == X3dnaBPType::cDUM) {
    return "c.-M";
  } else if (type == X3dnaBPType::cmPM) {
    return "cm+M";
  } else if (type == X3dnaBPType::cmUM) {
    return "cm-M";
  } else if (type == X3dnaBPType::DDDD) {
    return "....";
  } else if (type == X3dnaBPType::cmUW) {
    return "cm-W";
  } else if (type == X3dnaBPType::tMUm) {
    return "tM-m";
  } else if (type == X3dnaBPType::cDUW) {
    return "c.-W";
  } else if (type == X3dnaBPType::cMPm) {
    return "cM+m";
  } else if (type == X3dnaBPType::cMUm) {
    return "cM-m";
  } else if (type == X3dnaBPType::cDDD) {
    return "c...";
  } else if (type == X3dnaBPType::tWPm) {
    return "tW+m";
  } else if (type == X3dnaBPType::cDPm) {
    return "c.+m";
  } else if (type == X3dnaBPType::tmPm) {
    return "tm+m";
  } else if (type == X3dnaBPType::tWPD) {
    return "tW+.";
  } else if (type == X3dnaBPType::tmPW) {
    return "tm+W";
  } else if (type == X3dnaBPType::tDDD) {
    return "t...";
  } else if (type == X3dnaBPType::cWUD) {
    return "cW-.";
  } else if (type == X3dnaBPType::cWUM) {
    return "cW-M";
  } else if (type == X3dnaBPType::tDUW) {
    return "t.-W";
  } else if (type == X3dnaBPType::tMPM) {
    return "tM+M";
  } else if (type == X3dnaBPType::tDUM) {
    return "t.-M";
  } else if (type == X3dnaBPType::cMUD) {
    return "cM-.";
  } else if (type == X3dnaBPType::cWUm) {
    return "cW-m";
  } else if (type == X3dnaBPType::tDPm) {
    return "t.+m";
  } else if (type == X3dnaBPType::tMUD) {
    return "tM-.";
  } else if (type == X3dnaBPType::cmPW) {
    return "cm+W";
  } else if (type == X3dnaBPType::cMPM) {
    return "cM+M";
  } else if (type == X3dnaBPType::cmPD) {
    return "cm+.";
  } else if (type == X3dnaBPType::cmUD) {
    return "cm-.";
  } else if (type == X3dnaBPType::cDUD) {
    return "c.-.";
  } else if (type == X3dnaBPType::cWPW) {
    return "cW+W";
  } else if (type == X3dnaBPType::tDUD) {
    return "t.-.";
  } else if (type == X3dnaBPType::tDPW) {
    return "t.+W";
  } else if (type == X3dnaBPType::tmUm) {
    return "tm-m";
  } else if (type == X3dnaBPType::cWPD) {
    return "cW+.";
  } else if (type == X3dnaBPType::tmPD) {
    return "tm+.";
  } else if (type == X3dnaBPType::tDPD) {
    return "t.+.";
  } else if (type == X3dnaBPType::cDPD) {
    return "c.+.";
  } else if (type == X3dnaBPType::tDUm) {
    return "t.-m";
  } else if (type == X3dnaBPType::tDPM) {
    return "t.+M";
  }
  // added by CJ
  else if (type == X3dnaBPType::tWUD) {
    return "tW-.";
  } else if (type == X3dnaBPType::tmUW) {
    return "tm-W";
  } else if (type == X3dnaBPType::tMUM) {
    return "tM-M";
  } else if (type == X3dnaBPType::tMPD) {
    return "tM+.";
  } else if (type == X3dnaBPType::cDPW) {
    return "c.+W";
  } else if (type == X3dnaBPType::tmPM) {
    return "tm+M";
  } else if (type == X3dnaBPType::tWUm) {
    return "tW-m";
  } else if (type == X3dnaBPType::cWPm) {
    return "cW+m";
  } else if (type == X3dnaBPType::tmUD) {
    return "tm-.";
  } else if (type == X3dnaBPType::tWPM) {
    return "tW+M";
  } else if (type == X3dnaBPType::DWPm) {
    return ".W+m";
  } else if (type == X3dnaBPType::tMPW) {
    return "tM+W";
  } else if (type == X3dnaBPType::DDPm) {
    return "..+m";
  } else if (type == X3dnaBPType::tWUW) {
    return "tW-W";
  } else if (type == X3dnaBPType::cmPm) {
    return "cm+m";
  } else if (type == X3dnaBPType::DWUm) {
    return ".W-m";
  } else if (type == X3dnaBPType::DMPm) {
    return ".M+m";
  } else if (type == X3dnaBPType::DWPM) {
    return ".W+M";
  } else if (type == X3dnaBPType::DMPM) {
    return ".M+M";
  } else if (type == X3dnaBPType::DmPW) {
    return ".m+W";
  } else if (type == X3dnaBPType::DWUM) {
    return ".W-M";
  } else if (type == X3dnaBPType::DmPm) {
    return ".m+m";
  } else if (type == X3dnaBPType::DDUM) {
    return "..-M";
  } else if (type == X3dnaBPType::DMUm) {
    return ".M-m";
  } else if (type == X3dnaBPType::DDUm) {
    return "..-m";
  } else if (type == X3dnaBPType::DMPW) {
    return ".M+W";
  } else if (type == X3dnaBPType::DMPD) {
    return ".M+.";
  } else if (type == X3dnaBPType::DMUM) {
    return ".M-M";
  } else if (type == X3dnaBPType::DmUm) {
    return ".m-m";
  } else if (type == X3dnaBPType::DMUW) {
    return ".M-W";
  } else if (type == X3dnaBPType::DWUD) {
    return ".W-.";
  }

  else {
    throw X3dnaException("unknown x3dna bp type");
  }
}

/*
X3dna::X3Motifs get_motifs(String const &pdb_path) {

  auto dssr_file_path = _get_dssr_file_path(pdb_path);
  auto sections = _divide_dssr_file_into_sections(dssr_file_path);
  X3dna::X3Motifs all_motifs;

  StringStringMap types;
  types["hairpin"] = "HAIRPIN";
  types["bulges"] = "TWOWAY";
  types["internal"] = "TWOWAY";
  types["junction"] = "NWAY";
  types["non-loop"] = "SSTRAND";

  for (auto const &kv : types) {
    if (sections.find(kv.first) != sections.end()) {
      auto motifs = _parse_dssr_section(sections[kv.first], kv.second);
      for (auto const &m : motifs) {
        all_motifs.push_back(m);
      }
    }
  }
  if (sections.find("stems") != sections.end()) {
    auto motifs = _parse_dssr_helix_section(sections["stems"]);
    for (auto const &m : motifs) {
      all_motifs.push_back(m);
    }
  }

  return all_motifs;
}
 */
/*
X3dna::X3Motifs X3dna::_parse_dssr_section(Strings const &section,
                                           String const &mtype) {
  X3Motifs motifs;
  X3Residues seen_res;
  int count = 0;
  for (auto const &l : section) {
    auto spl = _split_over_white_space(l);
    if (spl.size() == 0) {
      continue;
    }
    try {
      if (spl[0].length() < 3 || spl[0].substr(0, 3) != "nts") {
        continue;
      }
    } catch (...) {
      continue;
    }
    if (spl.size() < 3) {
      continue;
    }

    auto res_strs = base::string::split(spl[2], ",");
    X3Residues res;
    for (auto const &res_str : res_strs) {
      auto res_obj = _parse_dssr_res_str(res_str);
      res.push_back(res_obj);
    }
    count = 0;
    for (auto const &r : res) {
      for (auto const &r2 : seen_res) {
        if (r == r2) {
          count += 1;
          break;
        }
      }
    }
    if (count == res.size()) {
      continue;
    }
    for (auto const &r : res) {
      seen_res.push_back(r);
    }
    motifs.push_back(X3Motif{res, mtype});
  }

  return motifs;
}
*/
/*
X3dna::X3Motifs X3dna::_parse_dssr_helix_section(Strings const &section) {

  X3Motifs motifs;
  X3Residues res;
  int i = 0;
  for (auto const &l : section) {
    auto spl = _split_over_white_space(l);
    if (spl.size() == 0) {
      continue;
    }
    try {
      i = std::stoi(spl[0]);
    } catch (...) {
      continue;
    }
    if (i == 1 && res.size() > 0) {
      motifs.push_back(X3Motif{res, "HELIX"});
      res = X3Residues();
    }
    res.push_back(_parse_dssr_res_str(spl[1]));
    res.push_back(_parse_dssr_res_str(spl[2]));
  }

  if (res.size() > 0) {
    motifs.push_back(X3Motif{res, "HELIX"});
  }

  return motifs;
}
*/

struct X3Basepair {
  // X3Residue res1, res2;
  math::Vector3 d;
  math::Matrix3x3 r;
  X3dnaBPType bp_type;
};

String X3dna::X3Basepair::to_string() const {
  auto ss = std::stringstream();
  ss << get_str_from_x3dna_type(bp_type) << "|" << res1.num << "|" << res2.num
     << "|" << d << "|" << math::Matrix3x3::matrix_to_str(r);
  return ss.str();
}

/*
void json_cleanup() {
  auto cleanup_cmd =
      String{base::x3dna_path()} + String{"/bin/x3dna-dssr --clean 2>"};
  std::system(cleanup_cmd.c_str());
}
*/

String compare_bps(X3dna::X3Basepairs &lhs, X3dna::X3Basepairs &rhs) {
  // first, we want to build the maps
  auto left_map = std::map<triplet<int, int, int>, X3dna::X3Basepair>();
  auto right_map = std::map<triplet<int, int, int>, X3dna::X3Basepair>();

  for (auto &bp : lhs) {
    auto key = triplet<int, int, int>();
    key.first = round(100. * bp.d.get_x());
    key.second = round(100. * bp.d.get_y());
    key.third = round(100. * bp.d.get_z());
    left_map[key] = bp;
  }
  for (auto &bp : rhs) {
    auto key = triplet<int, int, int>();
    key.first = round(100. * bp.d.get_x());
    key.second = round(100. * bp.d.get_y());
    key.third = round(100. * bp.d.get_z());
    right_map[key] = bp;
  }
  // quick check that there are no basepair repeats
  if (left_map.size() != lhs.size()) {
    throw std::runtime_error("error, redudant residue numbers in lhs bps");
  }

  if (right_map.size() != rhs.size()) {
    throw std::runtime_error("error, redudant residue numbers in rhs bps");
  }
  auto matches(0);
  auto lhs_it = left_map.begin();
  const auto lhs_end = left_map.end();
  while (lhs_it != lhs_end) {
    // check if right_map has the iterator
    auto rhs_it = right_map.find(lhs_it->first);
    // if so, then check if the two are equivalent. delete if so
    if (rhs_it != right_map.end()) {
      // need to check that the reference frame, origin and bp_type are
      // identical
      const auto same_origin_flag =
          roughly_equal(rhs_it->second.d, lhs_it->second.d, 0.1);
      const auto same_frame_flag =
          roughly_equal(rhs_it->second.r, lhs_it->second.r, 0.1) ||
          roughly_equal(rhs_it->second.r.get_flip_orientation(),
                        lhs_it->second.r, 0.1);
      const auto same_bp_type =
          rhs_it->second.bp_type == lhs_it->second.bp_type;
      if (same_origin_flag && same_frame_flag && same_bp_type) {
        lhs_it = left_map.erase(lhs_it);
        right_map.erase(rhs_it);
        ++matches;
      } else {
        // otherwise iterate
        ++lhs_it;
      }
    } else {
      ++lhs_it;
    }
  }
  auto summary = std::to_string(matches) + "," +
                 std::to_string(left_map.size()) + "," +
                 std::to_string(right_map.size()) + ",";
  // check the sizes of the maps... anything leftotver means they are not unique
  // to that method
  if (!left_map.empty()) {
    auto diffs = Strings{};
    for (auto &kv : left_map) {
      diffs.push_back(kv.second.to_string());
    }
    summary += base::string::join(diffs, "_");
  }
  summary += ",";
  if (!right_map.empty()) {
    auto diffs = Strings{};
    for (auto &kv : right_map) {
      diffs.push_back(kv.second.to_string());
    }
    summary += base::string::join(diffs, "_");
  }
  summary += "\n";

  return summary;
}
} // namespace x3dna
} // namespace util

/// @brief - ana_fncs.cpp

enum {
  TOR_XXX = 0,
  TOR_ALPHA,
  TOR_BETA,
  TOR_GAMMA,
  TOR_DELTA,
  TOR_EPSILON,
  TOR_ZETA,
  TOR_CHI,
  TOR_ETA,
  TOR_THETA,
  TOR_ETA1,
  TOR_THETA1,
  TOR_ETA2,
  TOR_THETA2,
  TOR_V0,
  TOR_V1,
  TOR_V2,
  TOR_V3,
  TOR_V4,
  TOR_TM,
  TOR_PHASE
};

enum {
  NT_NUM = 0,
  NT_P,
  NT_O5,
  NT_C5,
  NT_C4,
  NT_C3,
  NT_O3,
  NT_C2,
  NT_C1,
  NT_O4,
  NT_N,
  NT_C,
  O3P_LKG
};

char *nt_atoms[] = {"XXXX", " P  ", " O5'", " C5'", " C4'", " C3'", " O3'",
                    " C2'", " C1'", " O4'", "N9N1", "C4C2", "OPLK"};

/* for calculating overlap area between polygons */
#define MNPOLY 1000 /* maximum size of a polygon */

typedef struct {
  double x;
  double y;
} point;

typedef struct {
  point ip;
  point rx;
  point ry;
  long inside;
} vertex;

static double get_chi360(double chi) { return (chi > 0) ? chi : chi + 360; }

static long in_trans(double chi360) { return dval_in_range(chi360, 165, 315); }

static void check_chi(double chi, char *fmt, char *bstr, char *ostr) {
  parcat(ostr, chi, fmt, bstr);

  /* see http://www.fli-leibniz.de/ImgLibDoc/nana/chi.gif */
  if (chi > EMPTY_CRITERION) {
    double chi360;
    chi360 = get_chi360(chi);

    if (dval_in_range(chi360, 45, 95)) /* [60..80], with 15 allowance */
      strcat(ostr, " syn  ");
    else if (in_trans(chi360)) /* [-60, 180], with 15 allowance: [165, 315] */
      strcat(ostr, " anti ");
    else
      strcat(ostr, "      ");
  } else
    strcat(ostr, " ---  ");
}

/* http://www.imb-jena.de/Piet/help/backbone.html
 *   BI: (epsilon-zeta)	= -160 ... +20
 *  BII: (epsilon-zeta)	=  +20 ... +200 */

static void check_BI_BII(double epsilon, double zeta, char *bstr, char *fmt,
                         char *ostr) {
  if (epsilon < EMPTY_CRITERION || zeta < EMPTY_CRITERION) {
    parcat(ostr, EMPTY_NUMBER, fmt, bstr);
    strcat(ostr, "  ---"); /* no data for classification */
  } else {
    double d;
    if (epsilon < 0)
      epsilon += 360;
    if (zeta < 0)
      zeta += 360;
    d = epsilon - zeta;
    parcat(ostr, d, fmt, bstr);

    if (d < 0)
      d += 360;

    strcat(ostr, dval_in_range(d, 20, 200) ? "  BII" : "  BI");
  }
}

static void output_header_sugar_torsion(FILE *fp) {
  fprintf(fp, "Sugar conformational parameters: \n\n"
              "Note: v0: C4'-O4'-C1'-C2'\n"
              "      v1: O4'-C1'-C2'-C3'\n"
              "      v2: C1'-C2'-C3'-C4'\n"
              "      v3: C2'-C3'-C4'-O4'\n"
              "      v4: C3'-C4'-O4'-C1'\n\n"
              "      tm: the amplitude of pucker\n"
              "      P:  the phase angle of pseudorotation\n\n");
}

static void output_header_ss_Zp_Dp(FILE *fp) {
  fprintf(
      fp,
      "    ssZp: single-stranded Zp, defined as the z-coordinate of the 3'\n"
      "            phosphorus atom (P) expressed in the standard reference\n"
      "            frame of preceding base; the value is POSITIVE when P lies\n"
      "            on the +z-axis side (base in anti conformation); NEGATIVE\n"
      "            if P is on the -z-axis side (base in syn conformation)\n"
      "      Dp: perpendicular distance of the 3' P atom to the glycosydic "
      "bond\n"
      "            [as per the MolProbity paper of Richardson et al. "
      "(2010)]\n\n");
}

static void add_sugar_pucker(double phase_angle, char *bstr, char *ostr) {
  static char *sugar_pucker[10] = {
      "C3'-endo", "C4'-exo ", "O4'-endo", "C1'-exo ", "C2'-endo",
      "C3'-exo ", "C4'-endo", "O4'-exo ", "C1'-endo", "C2'-exo "};

  char temp[BUF512];
  long m;

  if (phase_angle > EMPTY_CRITERION) {
    m = (long)floor(phase_angle / 36.0);
    sprintf(temp, "%12s", sugar_pucker[m]);
    strcat(ostr, temp);
  } else
    strcat(ostr, bstr);
}

static void output_header_bb_torsion(FILE *fp) {
  fprintf(fp, "Main chain and chi torsion angles: \n\n"
              "Note: alpha:   O3'(i-1)-P-O5'-C5'\n"
              "      beta:    P-O5'-C5'-C4'\n"
              "      gamma:   O5'-C5'-C4'-C3'\n"
              "      delta:   C5'-C4'-C3'-O3'\n"
              "      epsilon: C4'-C3'-O3'-P(i+1)\n"
              "      zeta:    C3'-O3'-P(i+1)-O5'(i+1)\n\n"
              "      chi for pyrimidines(Y): O4'-C1'-N1-C2\n"
              "          chi for purines(R): O4'-C1'-N9-C4\n\n");
}

static void output_header_BI_BII_chi_syn_anti(FILE *fp) {
  fprintf(fp, "          chi in [165, -45(315)] for anti conformation\n"
              "                 chi in [45, 95] for syn conformation\n\n"
              "          e-z: epsilon - zeta\n"
              "              BI:  e-z = [-160, +20]\n"
              "              BII: e-z = [+20, +200]\n\n");
}

static double get_ss_Zp(double *P_xyz, double *o_xyz, double *n_xyz) {
  double dd[4];
  ddxyz(o_xyz, P_xyz, dd);
  return dot(dd, n_xyz);
}

static double idx2torsion(long *idx, double **xyz, long chk_lkg) {
  double d = EMPTY_NUMBER, **xyz4;
  long i;
  xyz4 = dmatrix(1, 4, 1, 3);

  for (i = 1; i <= 4; i++) {
    if (!idx[i])
      break;
    cpxyz(xyz[idx[i]], xyz4[i]);
  }
  if (i > 4)
    d = chk_lkg ? torsion(xyz4) : torsion2(xyz4);

  free_dmatrix(xyz4, 1, 4, 1, 3);

  return d;
}

void populate_nt_info(long num_residue, long **seidx, char **ResName,
                      char *ChainID, long *ResSeq, char **Miscs, char *bseq,
                      char **nt_info) {
  char idmsg[BUF512];
  long i, ib;

  for (i = 1; i <= num_residue; i++) {
    ib = seidx[i][1];
    base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], bseq[i], 1,
             idmsg);
    if (str_pmatch(idmsg, "....>"))
      strcpy(nt_info[i], idmsg + 5);
    else
      strcpy(nt_info[i], idmsg);
  }
}

void populate_nt_list(long num_residue, long **seidx, long *RY, char *bseq,
                      char **AtomName, double **xyz, long **nt_list) {
  long i, j, k, ib, ie;

  for (i = 1; i <= num_residue; i++) {
    if (RY[i] < 0)
      continue;
    ib = seidx[i][1];
    ie = seidx[i][2];

    for (j = 1; j <= NT_O4; j++)
      nt_list[i][j] = find_1st_atom(nt_atoms[j], AtomName, ib, ie, "");

    if (RY[i] == 1) { /* R: puRine */
      nt_list[i][NT_N] = find_1st_atom(" N9 ", AtomName, ib, ie, "");
      nt_list[i][NT_C] = find_1st_atom(" C4 ", AtomName, ib, ie, "");
    } else if (RY[i] == 0) {                  /* Y: pyRimidine */
      if (bseq[i] == 'P' || bseq[i] == 'p') { /* pseudo-U */
        nt_list[i][NT_N] = find_1st_atom(" C5 ", AtomName, ib, ie, "");
        nt_list[i][NT_C] = find_1st_atom(" C4 ", AtomName, ib, ie, "");
      } else {
        nt_list[i][NT_N] = find_1st_atom(" N1 ", AtomName, ib, ie, "");
        nt_list[i][NT_C] = find_1st_atom(" C2 ", AtomName, ib, ie, "");
      }
    }

    k = 0;
    for (j = 1; j <= NT_C; j++)
      if (nt_list[i][j])
        k++;
    nt_list[i][NT_NUM] = k;
  }

  for (i = 1; i < num_residue; i++) {
    ib = nt_list[i][NT_O3];
    ie = nt_list[i + 1][NT_P];
    if (ib && ie && within_limits(xyz[ib], xyz[ie], 0.8, O3P_UPPER))
      nt_list[i][O3P_LKG] = true;
  }
}

/* Get PDB data file name, output file name, and pairing information
   pair_num matrix is allocated here but de-allocated elsewhere */
long **read_input(char *inpfile, char *pdbfile, char *outfile, long *ds,
                  long *num_bp, long *ip, long *hetatm) {
  char str[BUF512];
  long i, j = 0, k, **pair_num;
  FILE *fp;

  fp = open_file(inpfile, "r");

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%s", pdbfile) != 1)
    fatal("cannot read pdbfile name: %s\n", pdbfile); /* PDB file name */

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%s", outfile) != 1)
    fatal("cannot read outfile name: %s\n", outfile); /* output file name */

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", ds) != 1)
    fatal("cannot read strand information\n"); /* ds is a pointer */
  if (*ds != 2 && *ds != 1)
    fatal("allowed options: 2--duplex and 1--single helix\n");

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", num_bp) != 1)
    fatal("cannot read number of base-pairs\n"); /* num_bp is a pointer */
  if (*num_bp <= 0)
    fatal("illegal number (%ld <= 0) of base-pairs\n", *num_bp);

  if (fgets(str, sizeof str, fp) == NULL || /* ip/hetatm are pointers */
      ((j = sscanf(str, "%ld %ld", ip, hetatm)) != 1 && j != 2)) {
    if (feof(fp)) {
      *ip = 0;
      *hetatm = 0;
    } else
      fatal("cannot read integer pair numbering/hetero indicator)\n");
  }

  if (j == 1)
    *hetatm = 0;
  pair_num = lmatrix(1, *ds + 1, 1, *num_bp);

  if (*ip) { /* user input */
    for (i = 1; i <= *num_bp; i++) {
      if (fgets(str, sizeof str, fp) == NULL)
        fatal("cannot read pair list\n");
      if (*ds == 2) { /* duplex */
        k = sscanf(str, "%ld %ld %ld", &pair_num[1][i], &pair_num[2][i],
                   &pair_num[3][i]);
        if (k != 2 && k != 3)
          fatal("two serial numbers required\n");
        if (k == 2 || (k == 3 && !pair_num[3][i] && pair_num[3][i] != 1 &&
                       pair_num[3][i] != 9))
          pair_num[3][i] = 0;
        if (pair_num[1][i] <= 0 || pair_num[2][i] <= 0)
          fatal("residue serial number %ld or %ld <= 0\n", pair_num[1][i],
                pair_num[2][i]);
      } else { /* single helix */
        if (sscanf(str, "%ld", &pair_num[1][i]) != 1)
          fatal("one serial number required: %s\n", str);
        if (pair_num[1][i] <= 0)
          fatal("residue serial number %ld <= 0\n", pair_num[1][i]);
      }
    }
  }

  close_file(fp);
  return pair_num;
}

void print_header(long ds, long num_bp, long num, char *pdbfile, FILE *fp) {
  time_t run_time;

  fprintf(fp, "    %s\n", Gvars.X3DNA_VER);
  print_sep(fp, '*', 76);
  fprintf(fp, "1. The list of the parameters given below correspond to"
              " the 5' to 3' direction\n   of strand I and 3' to 5' direction"
              " of strand II.\n\n");
  fprintf(fp, "2. All angular parameters, except for the phase angle"
              " of sugar pseudo-\n   rotation, are measured in degrees in"
              " the range of [-180, +180], and all\n"
              "   displacements are measured in Angstrom units.\n");

  print_sep(fp, '*', 76);
  fprintf(fp, "File name: %s\n", pdbfile);

  run_time = time(NULL);
  fprintf(fp, "Date and time: %s\n", ctime(&run_time));

  fprintf(fp, "Number of %s%s: %ld\n", (ds == 2) ? "base-pair" : "base",
          (num_bp == 1) ? "" : "s", num_bp);
  fprintf(fp, "Number of atoms: %ld\n", num);

  print_sep(fp, '*', 76);
  print_pdb_title(pdbfile, "*", fp);
}

static void output_atom_xyz(FILE *fp, long idx, double **xyz, char *bstr,
                            char *fmt) {
  long i;

  if (idx)
    for (i = 1; i <= 3; i++)
      fprintf(fp, fmt, xyz[idx][i]);
  else
    for (i = 1; i <= 3; i++)
      fprintf(fp, "%s", bstr);
}

void output_Borg_P_C1_C4(long num_residue, double **org, double **xyz,
                         long **nt_list, char **nt_info) {
  char *bstr = "    ---- ", *fmt = "%9.3f";
  long i, j;
  FILE *fp;

  fp = open_file("Borg_P_C1_C4.dat", "w");
  print_sep(fp, '*', 76);
  fprintf(fp,
          "xyz coordinates of the base origin, the P, C1' and C4' atoms\n\n");
  fprintf(
      fp,
      "              base      Ox       Oy       Oz       Px       Py       Pz"
      "      C1x      C1y      C1z      C4x      C4y      C4z\n");

  for (i = 1; i <= num_residue; i++) {
    if (nt_list[i][NT_NUM] < 3)
      continue;

    fprintf(fp, "%4ld %s", i, nt_info[i]);
    for (j = 1; j <= 3; j++)
      fprintf(fp, fmt, org[i][j]);
    output_atom_xyz(fp, nt_list[i][NT_P], xyz, bstr, fmt);
    output_atom_xyz(fp, nt_list[i][NT_C1], xyz, bstr, fmt);
    output_atom_xyz(fp, nt_list[i][NT_C4], xyz, bstr, fmt);
    fprintf(fp, "\n");
  }

  close_file(fp);
}

/* Indexes for sugar, chi torsion, C6-C8 etc atoms */
void atom_list(long ds, long num_bp, long **pair_num, long **seidx, long *RY,
               char **bp_seq, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, long **phos, long **c6_c8,
               long **sugar, long **chi) {
  char c2c4[5], c6c8[5], idmsg[BUF512], n1n9[5];
  long i, ib, ie, ioffset, j, c1, o4, rnum;

  for (i = 1; i <= ds; i++) {
    for (j = 1; j <= num_bp; j++) {
      rnum = pair_num[i][j];
      ib = seidx[rnum][1];
      ie = seidx[rnum][2];

      if (RY[rnum] < 0)
        fatal("Non-base residue: %s\n", ResName[ib]);
      get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
      /* backbone: P-O5'-C5'-C4'-C3'-O3' */
      phos[i][j] = find_1st_atom(" P  ", AtomName, ib, ie, idmsg);
      phos[i + 2][j] = find_1st_atom(" O1P", AtomName, ib, ie, idmsg);
      phos[i + 4][j] = find_1st_atom(" O2P", AtomName, ib, ie, idmsg);
      /* sugar: C4'-O4'-C1'-C2'-C3' */
      ioffset = (j - 1) * 5;
      sugar[i][ioffset + 1] = find_1st_atom(" C4'", AtomName, ib, ie, idmsg);
      sugar[i][ioffset + 5] = find_1st_atom(" C3'", AtomName, ib, ie, idmsg);
      o4 = find_1st_atom(" O4'", AtomName, ib, ie, idmsg);
      c1 = find_1st_atom(" C1'", AtomName, ib, ie, idmsg);
      sugar[i][ioffset + 2] = o4;
      sugar[i][ioffset + 3] = c1;
      sugar[i][ioffset + 4] = find_1st_atom(" C2'", AtomName, ib, ie, idmsg);

      /* chi(R): O4'-C1'-N9-C4; chi(Y): O4'-C1'-N1-C2 */
      ioffset = (j - 1) * 4;
      chi[i][ioffset + 1] = o4;
      chi[i][ioffset + 2] = c1;
      if (RY[rnum] == 1) {
        strcpy(n1n9, " N9 ");
        strcpy(c2c4, " C4 ");
        strcpy(c6c8, " C8 ");
      } else if (RY[rnum] == 0) {
        strcpy(n1n9, " N1 ");
        strcpy(c2c4, " C2 ");
        strcpy(c6c8, " C6 ");
        if (bp_seq[i][j] == 'P' || bp_seq[i][j] == 'p') {
          strcpy(n1n9, " C5 ");
          strcpy(c2c4, " C4 ");
        }
      }
      chi[i][ioffset + 3] = find_1st_atom(n1n9, AtomName, ib, ie, idmsg);
      chi[i][ioffset + 4] = find_1st_atom(c2c4, AtomName, ib, ie, idmsg);
      c6_c8[i][j] = find_1st_atom(c6c8, AtomName, ib, ie, idmsg);
    }
  }
}

static void get_mc_6_torsions(long num_residue, long *RY, long **bb,
                              double **xyz, double **nt_bb_torsion) {
  long i, idx[5];

  for (i = 1; i <= num_residue; i++) {
    if (RY[i] < 0) /* not a nt */
      continue;

    /* alpha: O3'(i - 1) - P(i) - O5'(i) - C5'(i) */
    idx[1] = (i > 1) ? bb[i - 1][6] : 0;
    idx[2] = bb[i][1];
    idx[3] = bb[i][2];
    idx[4] = bb[i][3];
    nt_bb_torsion[i][1] = idx2torsion(idx, xyz, true);

    /* beta: P(i) - O5'(i) - C5'(i) - C4'(i) */
    idx[1] = bb[i][1];
    idx[2] = bb[i][2];
    idx[3] = bb[i][3];
    idx[4] = bb[i][4];
    nt_bb_torsion[i][2] = idx2torsion(idx, xyz, true);
    /* gamma: O5'(i) - C5'(i) - C4'(i) - C3'(i) */
    idx[1] = bb[i][2];
    idx[2] = bb[i][3];
    idx[3] = bb[i][4];
    idx[4] = bb[i][5];
    nt_bb_torsion[i][3] = idx2torsion(idx, xyz, true);

    /* delta: C5'(i) - C4'(i) - C3'(i) - O3'(i) */
    idx[1] = bb[i][3];
    idx[2] = bb[i][4];
    idx[3] = bb[i][5];
    idx[4] = bb[i][6];
    nt_bb_torsion[i][4] = idx2torsion(idx, xyz, true);

    /* epsilon: C4'(i) - C3'(i) - O3'(i) - P(i + 1) */
    idx[1] = bb[i][4];
    idx[2] = bb[i][5];
    idx[3] = bb[i][6];
    idx[4] = (i < num_residue) ? bb[i + 1][1] : 0;
    nt_bb_torsion[i][5] = idx2torsion(idx, xyz, true);

    /* zeta: C3'(i) - O3'(i) - P(i + 1) - O5'(i + 1) */
    idx[1] = bb[i][5];
    idx[2] = bb[i][6];
    idx[3] = (i < num_residue) ? bb[i + 1][1] : 0;
    idx[4] = (i < num_residue) ? bb[i + 1][2] : 0;
    nt_bb_torsion[i][6] = idx2torsion(idx, xyz, true);
  }
}

void get_nt_torsion(long num_residue, double **org, double **xyz,
                    long **nt_list, double **nt_torsion) {
  static double Pconst;
  double phase_angle, **temp_xyz;
  long i, j, k, im1, ip1, idx[BUF32];

  Pconst = sin(PI / 5) + sin(PI / 2.5);
  temp_xyz = dmatrix(1, BUF32, 1, 3);

  init_dmatrix(nt_torsion, 1, num_residue, 1, BUF32, EMPTY_NUMBER);

  for (i = 1; i <= num_residue; i++) {
    im1 = i - 1;
    ip1 = i + 1;
    /* alpha: O3'(i-1) - P(i) - O5'(i) - C5'(i) */
    idx[1] = (i > 1) ? nt_list[im1][NT_O3] : 0;
    idx[2] = nt_list[i][NT_P];
    idx[3] = nt_list[i][NT_O5];
    idx[4] = nt_list[i][NT_C5];
    idx[5] = nt_list[i][NT_C4]; /* beta: P(i) - O5'(i) - C5'(i) - C4'(i) */
    idx[6] = nt_list[i][NT_C3]; /* gamma: O5'(i) - C5'(i) - C4'(i) - C3'(i) */
    idx[7] = nt_list[i][NT_O3]; /* delta: C5'(i) - C4'(i) - C3'(i) - O3'(i) */
    idx[8] = (i < num_residue)
                 ? nt_list[ip1][NT_P]
                 : 0; /* epsilon: C4'(i) - C3'(i) - O3'(i) - P(i+1) */
    idx[9] = (i < num_residue)
                 ? nt_list[ip1][NT_O5]
                 : 0; /* zeta: C3'(i) - O3'(i) - P(i+1) - O5'(i+1) */

    /* no need to check for O3P_LKG: implicit with true option below */
    nt_torsion[i][TOR_ALPHA] = idx2torsion(idx, xyz, true);
    nt_torsion[i][TOR_BETA] = idx2torsion(idx + 1, xyz, true);
    nt_torsion[i][TOR_GAMMA] = idx2torsion(idx + 2, xyz, true);
    nt_torsion[i][TOR_DELTA] = idx2torsion(idx + 3, xyz, true);
    nt_torsion[i][TOR_EPSILON] = idx2torsion(idx + 4, xyz, true);
    nt_torsion[i][TOR_ZETA] = idx2torsion(idx + 5, xyz, true);

    /* chi: O4'-C1'-N9-C4 (R); O4'-C1'-N1-C2 (Y) */
    idx[1] = nt_list[i][NT_O4];
    idx[2] = nt_list[i][NT_C1];
    idx[3] = nt_list[i][NT_N];
    idx[4] = nt_list[i][NT_C];
    nt_torsion[i][TOR_CHI] = idx2torsion(idx, xyz, true);

    /* eta:    C4'(i-1)-P(i)-C4'(i)-P(i+1)    --- C4'
       theta:  P(i)-C4'(i)-P(i+1)-C4'(i+1)
       eta':   C1'(i-1)-P(i)-C1'(i)-P(i+1)    --- C1'
       theta': P(i)-C1'(i)-P(i+1)-C1'(i+1)
       eta":   O(i-1)-P(i)-O(i)-P(i+1)        --- Origin of base
       theta": P(i)-O(i)-P(i+1)-O(i+1)  */
    if (i < num_residue && nt_list[i][O3P_LKG]) {
      idx[1] = (i > 1 && nt_list[im1][O3P_LKG]) ? nt_list[im1][NT_C4] : 0;
      idx[2] = nt_list[i][NT_P];
      idx[3] = nt_list[i][NT_C4];
      idx[4] = nt_list[ip1][NT_P];
      idx[5] = nt_list[ip1][NT_C4];
      nt_torsion[i][TOR_ETA] = idx2torsion(idx, xyz, false);
      nt_torsion[i][TOR_THETA] = idx2torsion(idx + 1, xyz, false);

      /* using C1' instead of C4' */
      idx[1] = (i > 1 && nt_list[im1][O3P_LKG]) ? nt_list[im1][NT_C1] : 0;
      idx[3] = nt_list[i][NT_C1];
      idx[5] = nt_list[ip1][NT_C1];
      nt_torsion[i][TOR_ETA1] = idx2torsion(idx, xyz, false);
      nt_torsion[i][TOR_THETA1] = idx2torsion(idx + 1, xyz, false);

      /* using origin of base reference frame */
      if (idx[2] && idx[4]) {
        if (i > 1 && nt_list[im1][O3P_LKG])
          copy_dvector(temp_xyz[1], org[im1], 1, 3);
        copy_dvector(temp_xyz[2], xyz[idx[2]], 1, 3);
        copy_dvector(temp_xyz[3], org[i], 1, 3);
        copy_dvector(temp_xyz[4], xyz[idx[4]], 1, 3);
        copy_dvector(temp_xyz[5], org[ip1], 1, 3);
        if (i > 1 && nt_list[im1][O3P_LKG])
          nt_torsion[i][TOR_ETA2] = torsion2(temp_xyz);
        nt_torsion[i][TOR_THETA2] = torsion2(temp_xyz + 1);
      }
    }

    /* sugar ring torsion angles v0 to v4 */
    idx[1] = nt_list[i][NT_C4];
    idx[2] = nt_list[i][NT_O4];
    idx[3] = nt_list[i][NT_C1];
    idx[4] = nt_list[i][NT_C2];
    idx[5] = nt_list[i][NT_C3];
    idx[6] = nt_list[i][NT_C4];
    idx[7] = nt_list[i][NT_O4];
    idx[8] = nt_list[i][NT_C1];
    k = 0;
    for (j = 0; j <= 4; j++) {
      nt_torsion[i][TOR_V0 + j] = idx2torsion(idx + j, xyz, true);
      if (idx[j + 1])
        k++;
    }
    if (k != 5)
      continue;

    /* phase angle and amplitude of pseudorotation */
    phase_angle = atan2(nt_torsion[i][TOR_V4] + nt_torsion[i][TOR_V1] -
                            nt_torsion[i][TOR_V3] - nt_torsion[i][TOR_V0],
                        2.0 * nt_torsion[i][TOR_V2] * Pconst);
    nt_torsion[i][TOR_TM] = nt_torsion[i][TOR_V2] / cos(phase_angle);

    phase_angle = rad2deg(phase_angle);
    if (phase_angle < 0)
      phase_angle += 360;
    nt_torsion[i][TOR_PHASE] = phase_angle;
  }

  free_dmatrix(temp_xyz, 1, BUF32, 1, 3);
}

void get_ss_Zp_Dp(long num_residue, double **org, double **orien, double **xyz,
                  long **nt_list, double **ss_Zp_Dp) {
  long i, ip1, idxP, idxN, idxC1;

  init_dmatrix(ss_Zp_Dp, 1, 2, 1, num_residue, EMPTY_NUMBER);

  for (i = 1; i <= num_residue; i++) {
    if (nt_list[i][NT_NUM] < 3)
      continue;

    ip1 = i + 1;
    if (i < num_residue && nt_list[i][O3P_LKG] && nt_list[ip1][NT_P]) {
      idxP = nt_list[ip1][NT_P];
      ss_Zp_Dp[1][i] = get_ss_Zp(xyz[idxP], org[i], orien[i] + 6);

      idxN = nt_list[i][NT_N];
      idxC1 = nt_list[i][NT_C1];
      if (idxC1 && idxN)
        ss_Zp_Dp[2][i] =
            get_point2line_perp_distance(xyz[idxP], xyz[idxC1], xyz[idxN]);
    }
  }
}

static long has_model_number(long num_residue, char **nt_info, long **nt_list) {
  long i;

  for (i = 1; i <= num_residue; i++) {
    if (nt_list[i][NT_NUM] < 3)
      continue;
    if (nt_info[i][4] == '>')
      return true;
  }

  return false;
}

void output_nt_torsion(long num_residue, char **nt_info, long **nt_list,
                       double **nt_torsion, double **ss_Zp_Dp, FILE *fp) {
  char *bstr = "    --- ", *fmt = "%8.1f";
  char ostr[BUF512];
  long i, j, with_model_number;

  with_model_number = has_model_number(num_residue, nt_info, nt_list);

  print_sep(fp, '*', 76);
  output_header_bb_torsion(fp);
  output_header_BI_BII_chi_syn_anti(fp);
  fprintf(fp,
          "              %sbase      chi A/S     alpha    beta   gamma   delta"
          "  epsilon   zeta     e-z BI/BII\n",
          with_model_number ? "     " : "");
  for (i = 1; i <= num_residue; i++) {
    if (nt_list[i][NT_NUM] < 3)
      continue;

    sprintf(ostr, "%4ld %s", i, nt_info[i]);
    check_chi(nt_torsion[i][TOR_CHI], fmt, bstr, ostr);

    for (j = TOR_ALPHA; j <= TOR_ZETA; j++)
      parcat(ostr, nt_torsion[i][j], fmt, bstr);
    check_BI_BII(nt_torsion[i][TOR_EPSILON], nt_torsion[i][TOR_ZETA], bstr, fmt,
                 ostr);
    fprintf(fp, "%s\n", ostr);
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "Pseudo (virtual) eta/theta torsion angles:\n\n");
  fprintf(fp, "Note: eta:    C4'(i-1)-P(i)-C4'(i)-P(i+1)\n"
              "      theta:  P(i)-C4'(i)-P(i+1)-C4'(i+1)\n\n"
              "      eta':   C1'(i-1)-P(i)-C1'(i)-P(i+1)\n"
              "      theta': P(i)-C1'(i)-P(i+1)-C1'(i+1)\n\n"
              "      eta\":   Borg(i-1)-P(i)-Borg(i)-P(i+1)\n"
              "      theta\": P(i)-Borg(i)-P(i+1)-Borg(i+1)\n\n");
  fprintf(fp,
          "              %sbase      eta   theta    eta'  theta'    eta\"  "
          "theta\"\n",
          with_model_number ? "     " : "");

  for (i = 1; i <= num_residue; i++) {
    if (nt_list[i][NT_NUM] < 3)
      continue;

    sprintf(ostr, "%4ld %s", i, nt_info[i]);
    for (j = TOR_ETA; j <= TOR_THETA2; j++)
      parcat(ostr, nt_torsion[i][j], fmt, bstr);
    fprintf(fp, "%s\n", ostr);
  }

  print_sep(fp, '*', 76);
  output_header_sugar_torsion(fp);
  output_header_ss_Zp_Dp(fp);
  fprintf(fp,
          "              %sbase       v0      v1      v2      v3      v4"
          "     tm       P    Puckering    ssZp     Dp\n",
          with_model_number ? "     " : "");
  for (i = 1; i <= num_residue; i++) {
    if (nt_list[i][NT_NUM] < 3)
      continue;

    sprintf(ostr, "%4ld %s", i, nt_info[i]);
    for (j = TOR_V0; j <= TOR_PHASE; j++)
      parcat(ostr, nt_torsion[i][j], fmt, bstr);
    add_sugar_pucker(nt_torsion[i][TOR_PHASE], bstr, ostr);

    if (nt_torsion[i][TOR_PHASE] < EMPTY_CRITERION)
      strcat(ostr, "    ");
    parcat(ostr, ss_Zp_Dp[1][i], "%8.2f", bstr);
    parcat(ostr, ss_Zp_Dp[2][i], "%8.2f", bstr);

    fprintf(fp, "%s\n", ostr);
  }
}

void get_nt_bb_torsion(double **nt_bb_torsion, long num_residue, long **seidx,
                       long *RY, char **AtomName, char **ResName, char *ChainID,
                       long *ResSeq, char **Miscs, double **xyz) {
  char idmsg[BUF512];
  char *bb_atoms[] = {" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " C1'"};
  long i, j, ib, ie, num_bb = sizeof(bb_atoms) / sizeof(bb_atoms[0]);
  long **bb;

  bb = lmatrix(1, num_residue, 1, num_bb);
  for (i = 1; i <= num_residue; i++) {
    if (RY[i] < 0)
      continue;

    ib = seidx[i][1];
    ie = seidx[i][2];
    get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
    for (j = 1; j <= num_bb; j++)
      bb[i][j] = find_1st_atom(bb_atoms[j - 1], AtomName, ib, ie, idmsg);
  }

  init_dmatrix(nt_bb_torsion, 1, num_residue, 1, 6, EMPTY_NUMBER);
  get_mc_6_torsions(num_residue, RY, bb, xyz, nt_bb_torsion);

  free_lmatrix(bb, 1, num_residue, 1, num_bb);
}

static void get_chi_torsions(long ds, long num_bp, long **chi, double **xyz,
                             double **chi_angle) {
  double **xyz4;
  long i, j, idx, ioffset, m;

  xyz4 = dmatrix(1, 4, 1, 3);

  for (i = 1; i <= ds; i++) {
    for (j = 1; j <= num_bp; j++) {
      ioffset = (j - 1) * 4;
      for (m = 1; m <= 4; m++) {
        idx = chi[i][ioffset + m];
        if (!idx)
          break;
        cpxyz(xyz[idx], xyz4[m]);
      }
      if (m > 4) /* all 4 indexes are okay  */
        chi_angle[i][j] = torsion(xyz4);
    }
  }

  free_dmatrix(xyz4, 1, 4, 1, 3);
}

static void get_sugar_torsions(long ds, long num_bp, long **sugar, double **xyz,
                               double **sugar_angle) {
  static double Pconst;
  static long vidx[5][4] = /* index of v0 to v4 */
      {{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 1}, {4, 5, 1, 2}, {5, 1, 2, 3}};
  double P_angle, **xyz4;
  long i, j, k, i5, idx, ioffset, m, o7;

  Pconst = sin(PI / 5) + sin(PI / 2.5);
  xyz4 = dmatrix(1, 4, 1, 3);

  for (i = 1; i <= ds; i++) {
    for (j = 1; j <= num_bp; j++) {
      /* sugar ring torsion angles v0 to v4 */
      i5 = 0;
      ioffset = (j - 1) * 5;
      o7 = (j - 1) * 7;
      for (k = 1; k <= 5; k++) {
        for (m = 1; m <= 4; m++) {
          idx = sugar[i][ioffset + vidx[k - 1][m - 1]];
          if (!idx)
            break;
          cpxyz(xyz[idx], xyz4[m]);
        }
        if (m > 4) { /* all 4 indexes are okay  */
          sugar_angle[i][o7 + k] = torsion(xyz4);
          i5++;
        }
      }

      /* phase angle and amplitude of pseudorotation */
      if (i5 == 5) {
        P_angle = atan2(sugar_angle[i][o7 + 5] + sugar_angle[i][o7 + 2] -
                            sugar_angle[i][o7 + 4] - sugar_angle[i][o7 + 1],
                        2 * sugar_angle[i][o7 + 3] * Pconst);
        sugar_angle[i][o7 + 6] = sugar_angle[i][o7 + 3] / cos(P_angle);
        P_angle = rad2deg(P_angle);
        if (P_angle < 0)
          P_angle = P_angle + 360;
        sugar_angle[i][o7 + 7] = P_angle;
      }
    }
  }

  free_dmatrix(xyz4, 1, 4, 1, 3);
}

static void print_duplex_torsions(long num_bp, long **pair_num, char **bp_seq,
                                  double **chi_angle, double **nt_bb_torsion,
                                  char *bstr, char *fmt, FILE *fp) {
  char str[BUF512];
  long i, j, m, idx, ds = 2;

  output_header_bb_torsion(fp);

  for (i = 1; i <= ds; i++) {
    (i == 1) ? fprintf(fp, "Strand I\n") : fprintf(fp, "Strand II\n");
    fprintf(fp,
            "  base    alpha    beta   gamma   delta  epsilon   zeta    chi\n");

    for (j = 1; j <= num_bp; j++) {
      sprintf(str, "%4ld %c ", j, bp_seq[i][j]);
      idx = pair_num[i][j]; /* residue number */
      for (m = 1; m <= 6; m++)
        parcat(str, nt_bb_torsion[idx][m], fmt, bstr);
      parcat(str, chi_angle[i][j], fmt, bstr);
      fprintf(fp, "%s\n", str);
    }

    if (i == 1)
      fprintf(fp, "\n");
  }
}

static void print_ss_torsions(long num_bp, long **pair_num, char **bp_seq,
                              double **chi_angle, double **nt_bb_torsion,
                              char *bstr, char *fmt, FILE *fp) {
  char str[BUF512];
  long i = 1, j, m, idx;

  output_header_bb_torsion(fp);

  fprintf(fp,
          "  base    alpha    beta   gamma   delta  epsilon   zeta    chi\n");

  for (j = 1; j <= num_bp; j++) {
    sprintf(str, "%4ld %c ", j, bp_seq[i][j]);
    idx = pair_num[i][j]; /* residue number */
    for (m = 1; m <= 6; m++)
      parcat(str, nt_bb_torsion[idx][m], fmt, bstr);
    parcat(str, chi_angle[i][j], fmt, bstr);
    fprintf(fp, "%s\n", str);
  }
}

static void print_duplex_sugar_conformation(long num_bp, char **bp_seq,
                                            double **sugar_angle, char *bstr,
                                            char *fmt, FILE *fp) {
  char str[BUF512];
  long i, j, k, ioffset, ds = 2;
  for (i = 1; i <= ds; i++) {
    (i == 1) ? fprintf(fp, "Strand I\n") : fprintf(fp, "Strand II\n");
    fprintf(fp, " base       v0      v1      v2      v3      v4"
                "      tm       P    Puckering\n");

    for (j = 1; j <= num_bp; j++) {
      sprintf(str, "%4ld %c ", j, bp_seq[i][j]);
      ioffset = (j - 1) * 7;
      for (k = 1; k <= 7; k++)
        parcat(str, sugar_angle[i][ioffset + k], fmt, bstr);
      add_sugar_pucker(sugar_angle[i][ioffset + 7], bstr, str);
      fprintf(fp, "%s\n", str);
    }

    if (i == 1)
      fprintf(fp, "\n");
  }
}

static void print_ss_sugar_conformation(long num_bp, char **bp_seq,
                                        double **sugar_angle, char *bstr,
                                        char *fmt, FILE *fp) {
  char str[BUF512];
  long i = 1, j, k, ioffset;

  fprintf(fp, " base       v0      v1      v2      v3      v4"
              "      tm       P    Puckering\n");

  for (j = 1; j <= num_bp; j++) {
    sprintf(str, "%4ld %c ", j, bp_seq[i][j]);

    ioffset = (j - 1) * 7;
    for (k = 1; k <= 7; k++)
      parcat(str, sugar_angle[i][ioffset + k], fmt, bstr);
    add_sugar_pucker(sugar_angle[i][ioffset + 7], bstr, str);

    fprintf(fp, "%s\n", str);
  }
}

void backbone_torsion(long ds, long num_bp, long **pair_num, char **bp_seq,
                      long **sugar, long **chi, double **xyz,
                      double **nt_bb_torsion, FILE *fp) {
  char *bstr = "    --- ", *fmt = "%8.1f";
  double **chi_angle, **sugar_angle;
  long num_bpx7;

  num_bpx7 = num_bp * 7;

  chi_angle = dmatrix(1, ds, 1, num_bp);
  sugar_angle = dmatrix(1, ds, 1, num_bpx7);

  /* initialize with EMPTY_NUMBER */
  init_dmatrix(chi_angle, 1, ds, 1, num_bp, EMPTY_NUMBER);
  init_dmatrix(sugar_angle, 1, ds, 1, num_bpx7, EMPTY_NUMBER);
  get_chi_torsions(ds, num_bp, chi, xyz, chi_angle);
  get_sugar_torsions(ds, num_bp, sugar, xyz, sugar_angle);

  print_sep(fp, '*', 76);
  if (ds == 2)
    print_duplex_torsions(num_bp, pair_num, bp_seq, chi_angle, nt_bb_torsion,
                          bstr, fmt, fp);
  else
    print_ss_torsions(num_bp, pair_num, bp_seq, chi_angle, nt_bb_torsion, bstr,
                      fmt, fp);

  print_sep(fp, '*', 76);
  output_header_sugar_torsion(fp);
  if (ds == 2)
    print_duplex_sugar_conformation(num_bp, bp_seq, sugar_angle, bstr, fmt, fp);
  else
    print_ss_sugar_conformation(num_bp, bp_seq, sugar_angle, bstr, fmt, fp);

  free_dmatrix(chi_angle, 1, ds, 1, num_bp);
  free_dmatrix(sugar_angle, 1, ds, 1, num_bpx7);
}

/* Calculate same strand P-P, C1'-C1' distances */
void p_c1_dist(long ds, long num_bp, char **bp_seq, long **phos, long **chi,
               double **xyz, long *bphlx, FILE *fp) {
  char *bstr = "       ---", *fmt = "%10.2f";
  char str[BUF512];
  double **c1_dist, **p_dist;
  long i, ia, ib, j, nbpm1;

  nbpm1 = num_bp - 1;

  p_dist = dmatrix(1, ds, 1, nbpm1);
  c1_dist = dmatrix(1, ds, 1, nbpm1);

  init_dmatrix(p_dist, 1, ds, 1, nbpm1, EMPTY_NUMBER);
  init_dmatrix(c1_dist, 1, ds, 1, nbpm1, EMPTY_NUMBER);

  for (i = 1; i <= ds; i++) {
    for (j = 1; j <= nbpm1; j++) {
      if (bphlx[j]) /* helix break */
        continue;
      ia = phos[i][j];
      ib = phos[i][j + 1];
      if (ia && ib)
        p_dist[i][j] = p1p2_dist(xyz[ia], xyz[ib]);
      ia = chi[i][(j - 1) * 4 + 2];
      ib = chi[i][j * 4 + 2];
      if (ia && ib)
        c1_dist[i][j] = p1p2_dist(xyz[ia], xyz[ib]);
    }
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "Same strand P--P and C1'--C1' virtual bond distances\n\n");

  fprintf(fp, "                 Strand I");
  if (ds == 2)
    fprintf(fp, "                          Strand II");
  fprintf(fp, "\n");

  fprintf(fp, "    step      P--P     C1'--C1'");
  if (ds == 2)
    fprintf(fp, "       step      P--P     C1'--C1'");
  fprintf(fp, "\n");

  for (i = 1; i <= nbpm1; i++) {
    j = 1; /* strand I */
    sprintf(str, "%4ld %c/%c", i, bp_seq[j][i], bp_seq[j][i + 1]);
    parcat(str, p_dist[j][i], fmt, bstr);
    parcat(str, c1_dist[j][i], fmt, bstr);
    fprintf(fp, "%s", str);

    if (ds == 2) {
      sprintf(str, "      %4ld %c/%c", i, bp_seq[ds][i], bp_seq[ds][i + 1]);
      parcat(str, p_dist[ds][i], fmt, bstr);
      parcat(str, c1_dist[ds][i], fmt, bstr);
      fprintf(fp, "%s", str);
    }
    fprintf(fp, "\n");
  }

  free_dmatrix(p_dist, 1, ds, 1, nbpm1);
  free_dmatrix(c1_dist, 1, ds, 1, nbpm1);
}

/* Get lambda angle and C1'-C1', C6-C8, N1-N9 distances */
void lambda_d3(long num_bp, char **bp_seq, long **chi, long **c6_c8,
               double **xyz, FILE *fp) {
  char str[BUF512], *bstr = "       ---", *fmt = "%10.1f";
  long i, j, ioffset, **c1_c1, **n1_n9;
  double vcc1[4], vcc2[4], vcn1[4], vcn2[4], **lambda_dist;

  c1_c1 = lmatrix(1, 2, 1, num_bp);
  n1_n9 = lmatrix(1, 2, 1, num_bp);
  lambda_dist = dmatrix(1, num_bp, 1, 5);

  init_dmatrix(lambda_dist, 1, num_bp, 1, 5, EMPTY_NUMBER);

  for (i = 1; i <= num_bp; i++) {
    ioffset = (i - 1) * 4;
    c1_c1[1][i] = chi[1][ioffset + 2];
    c1_c1[2][i] = chi[2][ioffset + 2];
    n1_n9[1][i] = chi[1][ioffset + 3];
    n1_n9[2][i] = chi[2][ioffset + 3];

    if (c1_c1[1][i] && c1_c1[2][i]) {
      ddxyz(xyz[c1_c1[2][i]], xyz[c1_c1[1][i]], vcc1);
      for (j = 1; j <= 3; j++)
        vcc2[j] = -vcc1[j];
      lambda_dist[i][3] = veclen(vcc1); /* C1'-C1' distance */
      if (n1_n9[1][i] && n1_n9[2][i]) {
        ddxyz(xyz[c1_c1[1][i]], xyz[n1_n9[1][i]], vcn1);
        ddxyz(xyz[c1_c1[2][i]], xyz[n1_n9[2][i]], vcn2);
        lambda_dist[i][1] = magang(vcc2, vcn1); /* lambda1 */
        lambda_dist[i][2] = magang(vcc1, vcn2); /* lambda2 */
      } else if (n1_n9[1][i] && !n1_n9[2][i]) {
        ddxyz(xyz[c1_c1[1][i]], xyz[n1_n9[1][i]], vcn1);
        lambda_dist[i][1] = magang(vcc2, vcn1); /* lambda1 */
      } else if (!n1_n9[1][i] && n1_n9[2][i]) {
        ddxyz(xyz[c1_c1[2][i]], xyz[n1_n9[2][i]], vcn2);
        lambda_dist[i][2] = magang(vcc1, vcn2); /* lambda2 */
      }
    }

    if (n1_n9[1][i] && n1_n9[2][i]) /* N1-N9 distance */
      lambda_dist[i][4] = p1p2_dist(xyz[n1_n9[1][i]], xyz[n1_n9[2][i]]);
    if (c6_c8[1][i] && c6_c8[2][i]) /* C6-C8 distance */
      lambda_dist[i][5] = p1p2_dist(xyz[c6_c8[1][i]], xyz[c6_c8[2][i]]);
  }
  print_sep(fp, '*', 76);
  fprintf(fp, "lambda: virtual angle between C1'-YN1 or C1'-RN9"
              " glycosidic bonds and the\n"
              "        base-pair C1'-C1' line\n\n"
              "C1'-C1': distance between C1' atoms for each base-pair\n"
              "RN9-YN1: distance between RN9-YN1 atoms for each base-pair\n"
              "RC8-YC6: distance between RC8-YC6 atoms for each base-pair\n");

  fprintf(fp,
          "\n    bp     lambda(I) lambda(II)  C1'-C1'   RN9-YN1   RC8-YC6\n");
  for (i = 1; i <= num_bp; i++) {
    sprintf(str, "%4ld %c%c%c", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
    for (j = 1; j <= 5; j++)
      parcat(str, lambda_dist[i][j], fmt, bstr);
    fprintf(fp, "%s\n", str);
  }

  /* write C1', RN9/YN1, RC8/YC6 xyz coordinates to "auxiliary.par" */
  print_axyz(num_bp, bp_seq, c1_c1, "C1'", xyz);
  print_axyz(num_bp, bp_seq, n1_n9, "RN9/YN1", xyz);
  print_axyz(num_bp, bp_seq, c6_c8, "RC8/YC6", xyz);

  free_lmatrix(c1_c1, 1, 2, 1, num_bp);
  free_lmatrix(n1_n9, 1, 2, 1, num_bp);
  free_dmatrix(lambda_dist, 1, num_bp, 1, 5);
}

/* Print xyz coordinates of P, C1', RN9/YN1 and RC8/YC6 */
void print_axyz(long num_bp, char **bp_seq, long **aidx, char *aname,
                double **xyz) {
  char *bstr = "    ---- ", *fmt = "%9.3f";
  long i, j;
  FILE *fc;

  fc = open_file(AUX_FILE, "a");
  print_sep(fc, '*', 76);
  fprintf(fc, "xyz coordinates of %s atoms\n\n", aname);

  fprintf(fc,
          "    bp        xI       yI       zI       xII      yII      zII\n");
  for (i = 1; i <= num_bp; i++) {
    fprintf(fc, "%4ld %c%c%c ", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
    if (aidx[1][i])
      for (j = 1; j <= 3; j++)
        fprintf(fc, fmt, xyz[aidx[1][i]][j]);
    else
      fprintf(fc, "%s%s%s", bstr, bstr, bstr);

    if (aidx[2][i])
      for (j = 1; j <= 3; j++)
        fprintf(fc, fmt, xyz[aidx[2][i]][j]);
    else
      fprintf(fc, "%s%s%s", bstr, bstr, bstr);

    fprintf(fc, "\n");
  }

  close_file(fc);
}

/* Calculate the correction angle for refined groove width */
static double gw_angle(long *idx, double *pvec12, double **xyz) {
  double pvec1[4], pvec2[4], pvecm[4];

  ddxyz(xyz[idx[2]], xyz[idx[1]], pvec1); /* strand I */
  ddxyz(xyz[idx[4]], xyz[idx[3]], pvec2); /* strand II */
  vec_norm(pvec1);
  vec_norm(pvec2);
  sumxyz(pvec1, pvec2, pvecm);

  return magang(pvecm, pvec12);
}

/* Groove width parameters based on El Hassan and Calladine (1998) */
void groove_width(long parallel, long num_bp, char **bp_seq, long **phos,
                  double **xyz, long *bphlx, FILE *fp) {
  char str[BUF512];
  char *bstr = "       ---", *fmt = "%10.1f"; /* groove width */
  char *bstr2 = "   ----", *fmt2 = "%7.1f";   /* P-P distance matrix */

  double anga, angb, dpa, dpb;
  double vecpa[4], vecpb[4];
  double **gwidth, **pdist;

  long iminor[5] = {BUF512, 1, -2, 2, -1};
  long imajor[3] = {BUF512, -2, 2};
  long idx[5], idxm1[5], idxp1[5], pp[5], ppa[5], ppb[5];
  long i, j, k, nbpm1, nset;
  long first_dist_num, items_per_line = 12, num_dist_sets;

  FILE *fc;

  nbpm1 = num_bp - 1;
  gwidth = dmatrix(1, nbpm1, 1, 4);
  init_dmatrix(gwidth, 1, nbpm1, 1, 4, EMPTY_NUMBER);

  /* method 1 is based on direct P-P distance
   *   minor: 0.5*[(P(i+1)-p(i-2))+(P(i+2)-p(i-1))]
   *   major: P(i-2)-p(i+2)
   * method 2 is method 1 plus a refinement
   *   minor: 0.5*[(P(i+1)-p(i-2))*sin(t1)+(P(i+2)-p(i-1))*sin(t2)]
   *   major: [P(i-2)-p(i+2)]*sin(t)
   */
  for (i = 1; i <= nbpm1; i++) {
    /* minor groove width */
    for (j = 1; j <= 4; j++) {
      idx[j] = i + iminor[j];
      if (!lval_in_range(idx[j], 1, nbpm1))
        break;
    }
    if (j > 4) {
      pp[1] = phos[1][idx[1] + 1];
      pp[3] = phos[1][idx[3] + 1];
      if (parallel) {
        pp[2] = phos[2][idx[2] + 1];
        pp[4] = phos[2][idx[4] + 1];
      } else {
        pp[2] = phos[2][idx[2]];
        pp[4] = phos[2][idx[4]];
      }

      if (pp[1] && pp[2] && pp[3] && pp[4]) {
        ddxyz(xyz[pp[2]], xyz[pp[1]], vecpa);
        ddxyz(xyz[pp[4]], xyz[pp[3]], vecpb);
        dpa = veclen(vecpa);
        dpb = veclen(vecpb);
        gwidth[i][1] = 0.5 * (dpa + dpb); /* method 1 */

        /* method 2 (refined P-P distance) */
        for (k = 1; k <= 4; k++) {
          idxm1[k] = idx[k] - 1;
          idxp1[k] = idx[k] + 1;
          if (idxm1[k] < 1 || idxp1[k] > nbpm1)
            break;
        }
        if (k > 4) {
          ppa[1] = phos[1][idxp1[1] + 1];
          ppa[2] = phos[1][idxm1[1] + 1];
          ppb[1] = phos[1][idxp1[3] + 1];
          ppb[2] = phos[1][idxm1[3] + 1];
          if (parallel) {
            ppa[3] = phos[2][idxp1[2] + 1];
            ppa[4] = phos[2][idxm1[2] + 1];
            ppb[3] = phos[2][idxp1[4] + 1];
            ppb[4] = phos[2][idxm1[4] + 1];
          } else {
            ppa[3] = phos[2][idxp1[2]];
            ppa[4] = phos[2][idxm1[2]];
            ppb[3] = phos[2][idxp1[4]];
            ppb[4] = phos[2][idxm1[4]];
          }
          if (ppa[1] && ppa[2] && ppa[3] && ppa[4] && ppb[1] && ppb[2] &&
              ppb[3] && ppb[4]) {
            anga = deg2rad(gw_angle(ppa, vecpa, xyz));
            angb = deg2rad(gw_angle(ppb, vecpb, xyz));
            gwidth[i][2] = 0.5 * (dpa * sin(anga) + dpb * sin(angb));
          }
        }
      }
    }
    /* major groove width */
    for (j = 1; j <= 2; j++) {
      idx[j] = i + imajor[j];
      if (!lval_in_range(idx[j], 1, nbpm1))
        break;
    }
    if (j > 2) {
      pp[1] = phos[1][idx[1] + 1];
      pp[2] = (parallel) ? phos[2][idx[2] + 1] : phos[2][idx[2]];
      if (pp[1] && pp[2]) {
        ddxyz(xyz[pp[2]], xyz[pp[1]], vecpa);
        dpa = veclen(vecpa);
        gwidth[i][3] = dpa; /* method 1 */
                            /* method 2 (refined P-P distance) */
        for (k = 1; k <= 2; k++) {
          idxm1[k] = idx[k] - 1;
          idxp1[k] = idx[k] + 1;
          if (idxm1[k] < 1 || idxp1[k] > nbpm1)
            break;
        }
        if (k > 2) {
          ppa[1] = phos[1][idxp1[1] + 1];
          ppa[2] = phos[1][idxm1[1] + 1];
          if (parallel) {
            ppa[3] = phos[2][idxp1[2] + 1];
            ppa[4] = phos[2][idxm1[2] + 1];
          } else {
            ppa[3] = phos[2][idxp1[2]];
            ppa[4] = phos[2][idxm1[2]];
          }
          if (ppa[1] && ppa[2] && ppa[3] && ppa[4]) {
            anga = deg2rad(gw_angle(ppa, vecpa, xyz));
            gwidth[i][4] = dpa * sin(anga);
          }
        }
      }
    }
  }
  print_sep(fp, '*', 76);
  fprintf(fp, "Minor and major groove widths: direct P-P distances "
              "and refined P-P distances\n   which take into account the "
              "directions of the sugar-phosphate backbones\n\n");
  fprintf(fp, "   (Subtract 5.8 Angstrom from the values to take account"
              " of the vdw radii\n    of the phosphate groups, and for"
              " comparison with FreeHelix and Curves.)\n\n");
  fprintf(fp, "Ref: M. A. El Hassan and C. R. Calladine (1998)."
              " ``Two Distinct Modes of\n     Protein-induced Bending"
              " in DNA.'' J. Mol. Biol., v282, pp331-343.\n\n");
  fprintf(fp, "                  Minor Groove        Major Groove\n"
              "                 P-P     Refined     P-P     Refined\n");

  /* taking into account helix breaks */
  for (i = 1; i <= nbpm1; i++) {
    if (bphlx[i]) {
      for (k = i - 3; k <= i + 3; k++) {
        if (!lval_in_range(k, 1, nbpm1))
          continue;
        for (j = 1; j <= 4; j++) /* for refined definition */
          if (!((k == i - 3 || k == i + 3) && (j == 1 || j == 3)))
            gwidth[k][j] = EMPTY_NUMBER;
      }
    }
  }
  for (i = 1; i <= nbpm1; i++) {
    sprintf(str, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
            bp_seq[2][i + 1], bp_seq[2][i]);
    for (j = 1; j <= 4; j++)
      parcat(str, gwidth[i][j], fmt, bstr);
    fprintf(fp, "%s\n", str);
  }
  /* P-P distance matrix */
  pdist = dmatrix(1, num_bp, 1, num_bp);
  init_dmatrix(pdist, 1, num_bp, 1, num_bp, EMPTY_NUMBER);

  for (i = 1; i <= num_bp; i++)
    for (j = 1; j <= num_bp; j++)
      if (phos[1][i] && phos[2][j])
        pdist[i][j] = p1p2_dist(xyz[phos[1][i]], xyz[phos[2][j]]);
  fc = open_file(AUX_FILE, "a");

  print_sep(fc, '*', 91);
  fprintf(fc, "Phosphorus-phosphorus distance in Angstroms\n\n");
  num_dist_sets = (long)ceil(num_bp / (double)items_per_line);
  for (nset = 1; nset <= num_dist_sets; nset++) {
    if (nset == num_dist_sets) {
      k = num_bp % items_per_line;
      if (!k)
        k = items_per_line;
    } else
      k = items_per_line;
    first_dist_num = (nset - 1) * items_per_line;
    fprintf(fc, "      ");
    for (i = first_dist_num + 1; i <= first_dist_num + k; i++)
      fprintf(fc, "%7ld", i);
    fprintf(fc, "\n");
    fprintf(fc, "         ");
    for (i = first_dist_num + 1; i <= first_dist_num + k; i++)
      fprintf(fc, "   %c   ", bp_seq[2][i]);
    fprintf(fc, "\n");

    for (i = 1; i <= num_bp; i++) {
      sprintf(str, "%4ld %c ", i, bp_seq[1][i]);
      for (j = first_dist_num + 1; j <= first_dist_num + k; j++)
        parcat(str, pdist[i][j], fmt2, bstr2);
      fprintf(fc, "%s\n", str);
    }
    if (nset != num_dist_sets)
      fprintf(fc, "\n");
  }
  close_file(fc);
  free_dmatrix(gwidth, 1, nbpm1, 1, 4);
  free_dmatrix(pdist, 1, num_bp, 1, num_bp);
}

void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch,
                          double opening) {
  static char *WC[9] = {WC_LIST};

  if (fabs(stretch) > 2.0 || fabs(opening) > 60)
    return;

  if (dval_in_range(fabs(shear), 1.8, 2.8))
    *bpid = 1; /* with WC geometry */
  if (fabs(shear) <= 1.8 && num_strmatch(bp, WC, 1, 8))
    *bpid = 2; /* WC */
}

/* Check if a base-pair is Watson-Crick:
   2: WC (1-below, plus |shear| <= 2.0, and base-pair sequence constraints)
   1: with correct WC geometry (i.e., x, y, z-axes in parallel directions)
   0: other cases, definitely non-WC (default) */
static void check_Watson_Crick(long num_bp, char **bp_seq, double **orien,
                               double **org, long *WC_info) {
  char bpi[3];
  double pars[7], o1[4], o2[4], morg[4], **r1, **r2, **mst;
  long i, j, k;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mst = dmatrix(1, 3, 1, 3);

  /* y- and z-axes of strand II base have been reversed */
  for (i = 1; i <= num_bp; i++) {
    j = (i - 1) * 9;
    sprintf(bpi, "%c%c", toupper((int)bp_seq[1][i]),
            toupper((int)bp_seq[2][i]));
    k = (bp_seq[0][i] ==
         '-') && /* anti-parallel: y- & z-axes already REVERSED */
        dot(&orien[1][j], &orien[2][j]) > 0.0 &&       /* x-axis */
        dot(&orien[1][j + 3], &orien[2][j + 3]) > 0.0; /* y-axis */
    if (k) {
      refs_right_left(i, orien, org, r1, o1, r2, o2);
      bpstep_par(r1, o1, r2, o2, pars, mst, morg);
      check_wc_wobble_pair(&WC_info[i], bpi, pars[1], pars[2], pars[6]);
    }
  }

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mst, 1, 3, 1, 3);
}

/* see also set_nmarkers() in 'find_pair.c' */
void set_chain_nmarkers019_to_symbols(long num, long *nmarkers,
                                      char *cmarkers) {
  char helix_begin = Gvars.CHAIN_MARKERS[0];
  char helix_middle = Gvars.CHAIN_MARKERS[1];
  char helix_end = Gvars.CHAIN_MARKERS[2];
  char isolated_bp = Gvars.CHAIN_MARKERS[3];
  long i, j, k, n, lastc;
  long *tcp, *temp, **seidx;

  lastc = strchr("0nofNOF", Gvars.CHAIN_MARKERS[4]) ? false : true;

  tcp = lvector(1, num);
  j = -1; /* index for isolated bps */
  k = 1;  /* index for for helices */
  for (i = 1; i <= num; i++) {
    if (nmarkers[i] == 0)
      tcp[i] = k;
    else if (nmarkers[i] == 9) {
      tcp[i] = k;
      k++;
    } else {
      tcp[i] = j;
      j--;
    }
  }

  temp = lvector(1, num);
  for (i = 1; i < num; i++)
    temp[i] = (tcp[i + 1] != tcp[i]) ? 1 : 0;
  temp[num] = 1;

  n = 0; /* get number of fragments */
  for (i = 1; i <= num; i++)
    if (temp[i])
      ++n;

  seidx = lmatrix(1, n, 0, 2); /* allocate spaces */
  n = 0;
  for (i = 1; i <= num; i++)
    if (temp[i])
      seidx[++n][2] = i;
  for (i = 2; i <= n; i++)
    seidx[i][1] = seidx[i - 1][2] + 1;
  seidx[1][1] = 1;

  for (i = 1; i <= n; i++)
    seidx[i][0] = seidx[i][2] - seidx[i][1] + 1;

  for (i = 1; i <= n; i++) {
    k = seidx[i][0];
    if (k == 1) { /* isolated bp */
      j = seidx[i][1];
      cmarkers[j] = isolated_bp;
    } else {
      for (j = seidx[i][1]; j <= seidx[i][2]; j++) {
        if (j == seidx[i][1])
          cmarkers[j] = helix_begin;
        else if (j == seidx[i][2])
          cmarkers[j] = helix_end;
        else
          cmarkers[j] = helix_middle;
      }
    }
  }

  if (!lastc && cmarkers[seidx[n][2]] != isolated_bp)
    cmarkers[seidx[n][2]] = helix_middle;
  free_lvector(tcp, 1, num);
  free_lvector(temp, 1, num);
  free_lmatrix(seidx, 1, n, 0, 2);
}

void get_bp_3char_symbols(long bp_type, char zdir, char *bp_sym) {
  sprintf(bp_sym, "%c%c%c", (bp_type == 2) ? '-' : '*',
          (bp_type > 0) ? '-' : '*', zdir);
}

/* Get the local reference frame for each base. Only the ring atoms are
   included in least-squares fitting */
void ref_frames(long ds, long num_bp, long **pair_num, char **bp_seq,
                long **seidx, long *RY, char **AtomName, char **ResName,
                char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                FILE *fp, double **orien, double **org, long *WC_info,
                long *str_type, long irna, long **o3p_brk) {
  static char *RingAtom[] = {RA_LIST};
  static char *rRingAtom[] =
      /* Babcock also uses C1* atom! */ {" C1'", RA_LIST};
  char bp_sym[BUF32], BDIR[BUF512], idmsg[BUF512], sidmsg[BUF512], spdb[BUF512];
  char *sChainID, **sAtomName, **sResName, **sMiscs, *cmarkers, *ss_markers;

  double orgi[4], vz[4];
  double **eRing_xyz, **fitted_xyz, **rms_fit, **sRing_xyz, **sxyz, **R;

  long i, ib, ie, ik, j, k, m, rnum, RingAtom_num, RA_NUM;
  long exp_katom, ioffset3, ioffset9, nmatch, snum, std_katom;
  long *sResSeq;

  irna ? get_BDIR(BDIR, "rAtomic_A.pdb") : get_BDIR(BDIR, "Atomic_A.pdb");

  RA_NUM = sizeof RingAtom / sizeof RingAtom[0];
  if (irna)
    RA_NUM++; /* One more C1' atom */
  rms_fit = dmatrix(1, ds, 1, num_bp);
  eRing_xyz = dmatrix(1, RA_NUM, 1, 3);
  sRing_xyz = dmatrix(1, RA_NUM, 1, 3);
  sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
  sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
  sChainID = cvector(1, NUM_RESIDUE_ATOMS);
  sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
  sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
  sMiscs = cmatrix(1, NUM_RESIDUE_ATOMS, 0, NMISC);
  fitted_xyz = dmatrix(1, RA_NUM, 1, 3);
  R = dmatrix(1, 3, 1, 3);

  for (i = 1; i <= ds; i++) {
    for (j = 1; j <= num_bp; j++) {
      rnum = pair_num[i][j];
      ib = seidx[rnum][1];
      ie = seidx[rnum][2];
      get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
      RingAtom_num = (RY[rnum] == 1) ? RA_NUM : RA_NUM - 3;
      set_std_base_pdb(BDIR, irna, bp_seq[i][j], spdb);
      snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                      sMiscs, 1, "*");
      sprintf(sidmsg, "in standard base: %s", spdb);

      nmatch = 0;
      for (k = 0; k < RingAtom_num; k++) {
        if (irna) {
          exp_katom = find_1st_atom(rRingAtom[k], AtomName, ib, ie, idmsg);
          std_katom = find_1st_atom(rRingAtom[k], sAtomName, 1, snum, sidmsg);
        } else {
          exp_katom = find_1st_atom(RingAtom[k], AtomName, ib, ie, idmsg);
          std_katom = find_1st_atom(RingAtom[k], sAtomName, 1, snum, sidmsg);
        }
        if (exp_katom && std_katom) {
          ++nmatch;
          cpxyz(xyz[exp_katom], eRing_xyz[nmatch]);
          cpxyz(sxyz[std_katom], sRing_xyz[nmatch]);
        }
      }
      rms_fit[i][j] =
          ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, orgi);
      ioffset9 = (j - 1) * 9;
      if (i == 2) {
        for (k = 1; k <= 3; k++) /* column-wise */
          vz[k] = R[k][3];
        if (dot(&orien[1][ioffset9 + 6], vz) < 0.0) {
          bp_seq[0][j] = '-'; /* anti-parallel */
          reverse_y_z_columns(R);
        } else
          bp_seq[0][j] = '+'; /* parallel */
      }
      cpxyz(orgi, org[i] + (j - 1) * 3);
      mst2orien(orien[i], ioffset9, R);
    }
  }

  if (ds == 2) {
    check_Watson_Crick(num_bp, bp_seq, orien, org, WC_info);
    ib = 0;
    ik = 0;
    if (num_bp == 1)
      init_dvector(orgi, 1, 3, 0.0);
    for (i = 1; i <= ds; i++) {
      for (j = 1; j <= num_bp; j++) {
        ioffset3 = (j - 1) * 3;
        if (j < num_bp)
          ddxyz(org[i] + ioffset3, org[i] + ioffset3 + 3, orgi);
        cpxyz(orien[i] + (j - 1) * 9 + 6, vz);
        if (!pair_num[3][j]) {
          ik++; /* non-breaks */
          if (dot(orgi, vz) < 0.0 && WC_info[j])
            ++ib; /* z-axis reversed */
        }
      }
    }

    /* most likely left-handed Z-DNA */
    if (ib && ib == ik) {
      *str_type = 1; /* with Z-axis reversed */
      for (i = 1; i <= ds; i++)
        for (j = 1; j <= num_bp; j++) {
          ioffset9 = (j - 1) * 9;
          negate_xyz(orien[i] + ioffset9);     /* reverse x-axis */
          negate_xyz(orien[i] + ioffset9 + 6); /* reverse z-axis */
        }
    }
    if (ib && ib != ik)
      *str_type = 2; /* unusual cases */
    if (ik != ds * num_bp)
      *str_type = *str_type + 10; /* more than one helix */
  }
  /* end of ds == 2 */
  /* write the least-squares fitting rms value */
  print_sep(fp, '*', 76);
  fprintf(fp, "RMSD of the bases");
  if (ds == 2)
    fprintf(fp, " (----- for WC bp, + for isolated bp, x for helix change)");

  fprintf(fp, "\n\n");
  fprintf(fp, "            Strand I");
  if (ds == 2)
    fprintf(fp, "                    Strand II          Helix");
  fprintf(fp, "\n");

  cmarkers = cvector(1, num_bp);
  ss_markers = cvector(1, num_bp);
  if (ds == 1)
    set_chain_nmarkers019_to_symbols(num_bp, o3p_brk[1], ss_markers);
  else
    set_chain_nmarkers019_to_symbols(num_bp, pair_num[ds + 1], cmarkers);

  for (i = 1; i <= num_bp; i++) {
    rnum = pair_num[1][i];
    ib = seidx[rnum][1];
    base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], bp_seq[1][i], 1,
             idmsg);
    if (ds == 1)
      fprintf(fp, "%4ld   (%5.3f) %s     %c\n", i, rms_fit[1][i], idmsg,
              ss_markers[i]);
    else {
      rnum = pair_num[2][i];
      ib = seidx[rnum][1];
      base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], bp_seq[2][i], 2,
               sidmsg);
      get_bp_3char_symbols(WC_info[i], bp_seq[0][i], bp_sym);
      fprintf(fp, "%4ld   (%5.3f) %s-%s-%s (%5.3f)     %c", i, rms_fit[1][i],
              idmsg, bp_sym, sidmsg, rms_fit[2][i], cmarkers[i]);
      if (Gvars.VERBOSE)
        fprintf(fp, " %c-%c", o3p_brk[1][i] ? 'x' : '-',
                o3p_brk[2][i] ? 'x' : '-');
      fprintf(fp, "\n");
    }
  }

  if (ds == 2) {
    k = 0;
    m = 0;
    for (i = 1; i <= num_bp; i++) {
      if (WC_info[i] != 2)
        k++;
      if (!WC_info[i])
        m++;
    }
    if (k)
      fprintf(fp,
              "\nNote: This structure contains %ld[%ld] non-Watson-Crick"
              " base-pair%s.\n",
              k, m, (k == 1) ? "" : "s");
  }
  free_dmatrix(rms_fit, 1, ds, 1, num_bp);
  free_dmatrix(eRing_xyz, 1, RA_NUM, 1, 3);
  free_dmatrix(sRing_xyz, 1, RA_NUM, 1, 3);
  free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq,
           sxyz, sMiscs);
  free_dmatrix(fitted_xyz, 1, RA_NUM, 1, 3);
  free_dmatrix(R, 1, 3, 1, 3);
  free_cvector(cmarkers, 1, num_bp);
  free_cvector(ss_markers, 1, num_bp);
}

/* Calculate step or base-pair parameters (using the CEHS scheme) */
void bpstep_par(double **rot1, double *org1, double **rot2, double *org2,
                double *pars, double **mst_orien, double *mst_org) {
  double phi, rolltilt;
  double hinge[4], mstx[4], msty[4], mstz[4], t1[4], t2[4];
  double **para_bp1, **para_bp2, **temp;
  long i, j;

  for (i = 1; i <= 3; i++) {
    t1[i] = rot1[i][3]; /* z1 */
    t2[i] = rot2[i][3]; /* z2 */
  }

  cross(t1, t2, hinge);
  rolltilt = magang(t1, t2);

  /* Handle the special cases where t1 and t2 are perfectly
     anti-parallel, thus hinge can not be defined by t1 and t2. The
     vector would be (0 0 0), corresponding to tilt^2 + roll^2 = 180
     Date: July 20, 2005 */
  if (veclen(hinge) < XEPS && (fabs(rolltilt - 180) < XEPS || rolltilt < XEPS))
    for (i = 1; i <= 3; i++)
      hinge[i] = rot1[i][1] + rot2[i][1] + rot1[i][2] + rot2[i][2];

  para_bp1 = dmatrix(1, 3, 1, 3);
  para_bp2 = dmatrix(1, 3, 1, 3);
  temp = dmatrix(1, 3, 1, 3);

  arb_rotation(hinge, -0.5 * rolltilt, temp);
  multi_matrix(temp, 3, 3, rot2, 3, 3, para_bp2);
  arb_rotation(hinge, 0.5 * rolltilt, temp);
  multi_matrix(temp, 3, 3, rot1, 3, 3, para_bp1);

  for (i = 1; i <= 3; i++) {
    mstz[i] = para_bp2[i][3]; /* also para_bp1(:,3) */
    t1[i] = para_bp1[i][2];   /* y1 */
    t2[i] = para_bp2[i][2];   /* y2 */
  }

  /* twist is the angle between the two y- or x-axes */
  pars[6] = vec_ang(t1, t2, mstz);

  /* get y- and x-axes */
  get_vector(t1, mstz, 0.5 * pars[6], msty);
  cross(msty, mstz, mstx);

  avexyz(org1, org2, mst_org);
  ddxyz(org1, org2, t1);
  x_y_z_2_mtx(mstx, msty, mstz, mst_orien);

  /* get the xyz displacement parameters */
  for (i = 1; i <= 3; i++) {
    pars[i] = 0.0;
    for (j = 1; j <= 3; j++)
      pars[i] += t1[j] * mst_orien[j][i];
  }

  /* phi angle is defined by hinge and msty */
  phi = deg2rad(vec_ang(hinge, msty, mstz));

  /* get roll and tilt angles */
  pars[5] = rolltilt * cos(phi);
  pars[4] = rolltilt * sin(phi);

  free_dmatrix(para_bp1, 1, 3, 1, 3);
  free_dmatrix(para_bp2, 1, 3, 1, 3);
  free_dmatrix(temp, 1, 3, 1, 3);
}

/* Calculate local helical parameters & its middle frame */
void helical_par(double **rot1, double *org1, double **rot2, double *org2,
                 double *pars, double **mst_orien, double *mst_org) {
  double AD_mag, phi, TipInc1, TipInc2, vlen;
  double axis_h[4], hinge1[4], hinge2[4], t1[4], t2[4];
  double AD_axis[4], org1_h[4], org2_h[4];
  double **rot1_h, **rot2_h, **temp;
  long i, j;

  for (i = 1; i <= 3; i++) {         /* column-wise */
    t1[i] = rot2[i][1] - rot1[i][1]; /* dx */
    t2[i] = rot2[i][2] - rot1[i][2]; /* dy */
  }
  cross(t1, t2, axis_h);
  vlen = veclen(axis_h);
  if (vlen < XEPS) { /* for twist = 0.0 */
    axis_h[1] = 0.0;
    axis_h[2] = 0.0;
    axis_h[3] = 1.0;
  } else
    for (i = 1; i <= 3; i++)
      axis_h[i] /= vlen;

  temp = dmatrix(1, 3, 1, 3);
  rot1_h = dmatrix(1, 3, 1, 3);
  for (i = 1; i <= 3; i++)
    t1[i] = rot1[i][3]; /* z1 */
  TipInc1 = magang(axis_h, t1);
  cross(axis_h, t1, hinge1);
  arb_rotation(hinge1, -TipInc1, temp);
  multi_matrix(temp, 3, 3, rot1, 3, 3, rot1_h);
  rot2_h = dmatrix(1, 3, 1, 3);
  for (i = 1; i <= 3; i++)
    t2[i] = rot2[i][3]; /* z2 */
  TipInc2 = magang(axis_h, t2);
  cross(axis_h, t2, hinge2);
  arb_rotation(hinge2, -TipInc2, temp);
  multi_matrix(temp, 3, 3, rot2, 3, 3, rot2_h);

  for (i = 1; i <= 3; i++) {
    t1[i] = rot1_h[i][1] + rot2_h[i][1]; /* x1 + x2 */
    t2[i] = rot1_h[i][2] + rot2_h[i][2]; /* y1 + y2 */
  }
  vec_norm(t1);
  vec_norm(t2);
  x_y_z_2_mtx(t1, t2, axis_h, mst_orien);

  for (i = 1; i <= 3; i++) {
    t1[i] = rot1_h[i][2]; /* y1_h */
    t2[i] = rot2_h[i][2]; /* y2_h */
  }
  pars[6] = vec_ang(t1, t2, axis_h);

  ddxyz(org1, org2, t2); /* org2-org1 */
  pars[3] = dot(t2, axis_h);

  phi = deg2rad(vec_ang(hinge1, t1, axis_h));
  pars[5] = TipInc1 * cos(phi);
  pars[4] = TipInc1 * sin(phi);

  for (i = 1; i <= 3; i++)
    t1[i] = t2[i] - pars[3] * axis_h[i];
  if (fabs(pars[6]) < HTWIST0) /* twist = 0.0: cf <xhelfunc> */
    for (i = 1; i <= 3; i++)
      org1_h[i] = org1[i] + 0.5 * t1[i];
  else {
    get_vector(t1, axis_h, 90 - pars[6] / 2, AD_axis);
    AD_mag = 0.5 * veclen(t1) / sin(deg2rad(pars[6] / 2));
    for (i = 1; i <= 3; i++)
      org1_h[i] = org1[i] + AD_mag * AD_axis[i];
  }

  for (i = 1; i <= 3; i++)
    org2_h[i] = org1_h[i] + pars[3] * axis_h[i];
  avexyz(org1_h, org2_h, mst_org);
  ddxyz(org1_h, org1, t1);

  for (i = 1; i <= 2; i++) {
    pars[i] = 0.0;
    for (j = 1; j <= 3; j++)
      pars[i] += t1[j] * rot1_h[j][i];
  }

  free_dmatrix(rot1_h, 1, 3, 1, 3);
  free_dmatrix(temp, 1, 3, 1, 3);
  free_dmatrix(rot2_h, 1, 3, 1, 3);
}

/* Print base-pair, step and helical parameters */
void print_par(char **bp_seq, long num_bp, long ich, long ishel, double **param,
               FILE *fp) {
  char *fmt = "%10.2f";
  double temp[7];
  long i, j;

  if (!lval_in_range(ich, 1, 4))
    fatal("wrong option [%ld] for printing parameters\n", ich);

  if (ich == 1) { /* base-pair parameters */
    fprintf(fp, "     bp        Shear    Stretch   Stagger"
                "    Buckle  Propeller  Opening\n");
    for (i = 1; i <= num_bp; i++) {
      fprintf(fp, " %4ld %c%c%c ", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
      for (j = 1; j <= 6; j++)
        fprintf(fp, fmt, param[i][j]);
      fprintf(fp, "\n");
    }

  } else {
    if (num_bp == 1)
      return;
    if (ishel)
      fprintf(fp, "    step       X-disp    Y-disp   h-Rise"
                  "     Incl.       Tip   h-Twist\n");
    else
      fprintf(fp, "    step       Shift     Slide      Rise"
                  "      Tilt      Roll     Twist\n");
    for (i = 1; i <= num_bp - 1; i++) {
      if (ich == 2) /* for base-pair step */
        fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
                bp_seq[2][i + 1], bp_seq[2][i]);
      else
        fprintf(fp, "%4ld  %c/%c ", i, bp_seq[ich - 2][i],
                bp_seq[ich - 2][i + 1]);
      for (j = 1; j <= 6; j++)
        fprintf(fp, fmt, param[i][j]);
      fprintf(fp, "\n");
    }
  }

  if (num_bp > 2) {
    j = (ich == 1) ? num_bp : num_bp - 1;
    fprintf(fp, "          ");
    print_sep(fp, '~', 60);
    fprintf(fp, "      ave.");
    ave_dmatrix(param, j, 6, temp);
    for (i = 1; i <= 6; i++)
      fprintf(fp, fmt, temp[i]);
    fprintf(fp, "\n");
    fprintf(fp, "      s.d.");
    std_dmatrix(param, j, 6, temp);
    for (i = 1; i <= 6; i++)
      fprintf(fp, fmt, temp[i]);
    fprintf(fp, "\n");
  }
}

/* Origin xyz coordinates followed by direction cosines of x-, y- & z-axes */
static void print_analyze_ref_frames(long ds, long num_bp, char **bp_seq,
                                     double *iorg, double *iorien,
                                     long **pair_num, char **nt_info) {
  FILE *fp;
  long i, j, ia, ib, ik, ioffset3, ioffset9;

  fp = open_file(REF_FILE, "w");

  if (ds == 1)
    fprintf(fp, "%5ld bases\n", num_bp);
  else
    fprintf(fp, "%5ld base-pairs\n", num_bp);

  for (i = 1; i <= num_bp; i++) {
    ia = pair_num[1][i];
    if (ds == 1)
      fprintf(fp, "... %5ld %c   # %s\n", i, bp_seq[1][i], nt_info[ia]);
    else {
      ib = pair_num[2][i];
      fprintf(fp, "... %5ld %c%c%c   # %s - %s\n", i, bp_seq[1][i],
              bp_seq[0][i], bp_seq[2][i], nt_info[ia], nt_info[ib]);
    }
    ioffset3 = (i - 1) * 3;
    fprintf(fp, "%10.4f %10.4f %10.4f  # origin\n", iorg[ioffset3 + 1],
            iorg[ioffset3 + 2], iorg[ioffset3 + 3]);
    ioffset9 = (i - 1) * 9;
    for (j = 1; j <= 3; j++) {
      ik = ioffset9 + (j - 1) * 3;
      fprintf(fp, "%10.4f %10.4f %10.4f  # %c-axis\n", iorien[ik + 1],
              iorien[ik + 2], iorien[ik + 3],
              (j == 1)   ? 'x'
              : (j == 2) ? 'y'
                         : 'z');
    }
  }

  close_file(fp);
}

static void single_helix(long num_bp, char **bp_seq, double **step_par,
                         double **heli_par, double **orien, double **org,
                         FILE *fp, long **pair_num, char **nt_info) {
  char str[BUF512];
  double **parmtx;
  long nbpm1;
  FILE *fstep, *fheli, *fchek;

  nbpm1 = num_bp - 1;

  parmtx = dmatrix(1, nbpm1, 1, 6);

  print_sep(fp, '*', 76);
  fprintf(fp, "Local base step parameters\n");
  parvec2mtx(step_par[1], nbpm1, parmtx);
  print_par(bp_seq, num_bp, 3, 0, parmtx, fp);

  /* step parameters for rebuilding */
  fstep = open_file(BPSTEP_FILE, "w");
  sprintf(str, "%4ld # bases\n%4ld # ***local step parameters***\n", num_bp,
          0L);
  strcat(str,
         "#      Shift     Slide     Rise      Tilt      Roll      Twist\n");
  print_ss_rebuild_pars(parmtx, num_bp, str, bp_seq, fstep);
  close_file(fstep);

  print_sep(fp, '*', 76);
  fprintf(fp, "Local base helical parameters\n");
  parvec2mtx(heli_par[1], nbpm1, parmtx);
  print_par(bp_seq, num_bp, 3, 1, parmtx, fp);

  /* helical parameters for rebuilding */
  fheli = open_file(HLXSTEP_FILE, "w");
  sprintf(str, "%4ld # bases\n%4ld # ***local helical parameters***\n", num_bp,
          1L);
  strcat(str,
         "#      X-disp    Y-disp    h-Rise    Incl.     Tip     h-Twist\n");
  print_ss_rebuild_pars(parmtx, num_bp, str, bp_seq, fheli);
  close_file(fheli);

  /* for checking out */
  fchek = open_file(AUX_FILE, "w");
  fprintf(fchek, "Reference frame: Origins (Ox, Oy, Oz) followed by the"
                 " direction cosines of the\n"
                 "                 X- (Xx, Yx, Zx), Y- (Yx, Yx, Yx),"
                 " and Z- (Zx, Zx, Zx) axes\n");
  print_sep(fchek, '*', 89);
  fprintf(fchek, "Local base reference frames\n");
  print_ref(bp_seq, num_bp, 2, org[1], orien[1], fchek);
  close_file(fchek);

  /* reference frame for reseting the structure */
  print_analyze_ref_frames(1, num_bp, bp_seq, org[1], orien[1], pair_num,
                           nt_info);

  free_dmatrix(parmtx, 1, num_bp, 1, 6);
}

void output_ave_std(long num, double **parcln, int dnum, char *fmt, FILE *fp) {
  double temp[7];
  int k = 10;
  long i;

  if (num <= 2)
    return;

  k += dnum;
  fprintf(fp, "%*s", k, " ");
  print_sep(fp, '~', 60);

  if (dnum)
    fprintf(fp, "%*s      ave.", dnum, " ");
  else
    fprintf(fp, "      ave.");
  ave_dmatrix(parcln, num, 6, temp);
  for (i = 1; i <= 6; i++)
    fprintf(fp, fmt, temp[i]);
  fprintf(fp, "\n");

  if (dnum)
    fprintf(fp, "%*s      s.d.", dnum, " ");
  else
    fprintf(fp, "      s.d.");
  std_dmatrix(parcln, num, 6, temp);
  for (i = 1; i <= 6; i++)
    fprintf(fp, fmt, temp[i]);
  fprintf(fp, "\n");
}

/* Print local base-pair step and helical parameter with helix breaks deleted */
void prt_stepstr(char **step_str, long num_step, long *bphlx, long ishel,
                 double **param, FILE *fp) {
  char *bstr = "      ----", *fmt = "%10.2f";
  double **parcln;
  long i, j, num = 0;

  if (!num_step) /* empty */
    return;

  parcln = dmatrix(1, num_step, 1, 6);
  if (ishel)
    fprintf(fp, "    step       X-disp    Y-disp   h-Rise"
                "     Incl.       Tip   h-Twist\n");
  else
    fprintf(fp, "    step       Shift     Slide      Rise"
                "      Tilt      Roll     Twist\n");
  for (i = 1; i <= num_step; i++) {
    fprintf(fp, "%4ld %s", i, step_str[i]);
    if (bphlx[i])
      for (j = 1; j <= 6; j++)
        fprintf(fp, "%s", bstr);
    else {
      num++;
      for (j = 1; j <= 6; j++) {
        fprintf(fp, fmt, param[i][j]);
        parcln[num][j] = param[i][j];
      }
    }
    fprintf(fp, "\n");
  }

  output_ave_std(num, parcln, 0, fmt, fp);

  free_dmatrix(parcln, 1, num_step, 1, 6);
}

/* Print local base-pair step and helical parameter with helix breaks
 * deleted. Some heavily kinked steps (e.g., pdt040) could be divided */
void prt_step_par(char **bp_seq, long num_bp, long *bphlx, long ishel,
                  double **param, FILE *fp) {
  char *bstr = "      ----", *fmt = "%10.2f";
  double temp[7], **parcln;
  long i, j, num = 0, nbpm1;

  if (num_bp == 1)
    return;

  nbpm1 = num_bp - 1;
  parcln = dmatrix(1, nbpm1, 1, 6);
  if (ishel)
    fprintf(fp, "    step       X-disp    Y-disp   h-Rise"
                "     Incl.       Tip   h-Twist\n");
  else
    fprintf(fp, "    step       Shift     Slide      Rise"
                "      Tilt      Roll     Twist\n");
  for (i = 1; i <= nbpm1; i++) {
    fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
            bp_seq[2][i + 1], bp_seq[2][i]);
    if (bphlx[i])
      for (j = 1; j <= 6; j++)
        fprintf(fp, "%s", bstr);
    else {
      num++;
      for (j = 1; j <= 6; j++) {
        fprintf(fp, fmt, param[i][j]);
        parcln[num][j] = param[i][j];
      }
    }
    fprintf(fp, "\n");
  }

  if (num > 2) {
    fprintf(fp, "          ");
    print_sep(fp, '~', 60);
    fprintf(fp, "      ave.");
    ave_dmatrix(parcln, num, 6, temp);
    for (i = 1; i <= 6; i++)
      fprintf(fp, fmt, temp[i]);
    fprintf(fp, "\n");
    fprintf(fp, "      s.d.");
    std_dmatrix(parcln, num, 6, temp);
    for (i = 1; i <= 6; i++)
      fprintf(fp, fmt, temp[i]);
    fprintf(fp, "\n");
  }
  free_dmatrix(parcln, 1, nbpm1, 1, 6);
}

static void double_helix(long num_bp, char **bp_seq, double **step_par,
                         double **heli_par, double **orien, double **org,
                         long *WC_info, FILE *fp, double **twist_rise,
                         double *mst_orien, double *mst_org, double *mst_orienH,
                         double *mst_orgH, long *bphlx, long istart, long istep,
                         long bz, long *str_type, long **pair_num,
                         char **nt_info) {
  char str[BUF512], *fmt = " %9.3f";
  char **step_str;
  double hfoi[4], mfoi[4], o1[4], o2[4];
  double *bp_org, *bp_orien, **bp_par, **bp_step_par, **bp_heli_par;
  double **mfi, **hfi, **r1, **r2;
  long num_step = 0, bz_junction = 0, z_step = 0;
  long i, ib, ie, ioffset3, ioffset9, j, nbpm1;
  FILE *fchek, *fheli, *fstep;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);
  hfi = dmatrix(1, 3, 1, 3);

  bp_par = dmatrix(1, num_bp, 1, 6);
  bp_org = dvector(1, num_bp * 3);
  bp_orien = dvector(1, num_bp * 9);

  print_sep(fp, '*', 76);
  fprintf(fp, "Origin (Ox, Oy, Oz) and mean normal vector"
              " (Nx, Ny, Nz) of each base-pair in\n"
              "   the coordinate system of the given structure\n\n");
  fprintf(
      fp,
      "      bp        Ox        Oy        Oz        Nx        Ny        Nz\n");

  for (i = 1; i <= num_bp; i++) {
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    bpstep_par(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

    fprintf(fp, " %4ld %c%c%c ", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
    for (j = 1; j <= 3; j++)
      fprintf(fp, fmt, mfoi[j]); /* origin */
    for (j = 1; j <= 3; j++)
      fprintf(fp, fmt, mfi[j][3]); /* base-pair normal */
    fprintf(fp, "\n");

    cpxyz(mfoi, bp_org + (i - 1) * 3);
    mst2orien(bp_orien, (i - 1) * 9, mfi);
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "Local base-pair parameters\n");
  print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

  bp_step_par = dmatrix(1, nbpm1, 1, 6);
  bp_heli_par = dmatrix(1, nbpm1, 1, 6);
  step_str = cmatrix(1, num_bp, 0, 5);
  ie = istart;
  for (i = istart; i <= nbpm1; i++) {
    if (istep > 0) { /* continuous 1 to 3; 2 to 4 etc */
      ib = i;
      ie = ib + istep;
    } else { /* next segment: 1 to 3, 3 to 5 etc */
      ib = ie;
      ie = ib - istep;
    }
    if (ie > num_bp)
      break;
    num_step++;
    sprintf(step_str[num_step], "%c%c/%c%c", bp_seq[1][ib], bp_seq[1][ie],
            bp_seq[2][ie], bp_seq[2][ib]);
    refs_i_j(ib, ie, bp_orien, bp_org, r1, o1, r2, o2);
    if (WC_info[i] == 2 && WC_info[i + 1] == 2) /* must be WC step */
      bz_check(r1, o1, r2, o2, bz, &bz_junction, &z_step);
    bpstep_par(r1, o1, r2, o2, bp_step_par[num_step], mfi, mfoi);
    helical_par(r1, o1, r2, o2, bp_heli_par[num_step], hfi, hfoi);

    /* this could be meaningless for non-dinucleotide steps */
    twist_rise[num_step][1] = bp_step_par[num_step][6];
    twist_rise[num_step][2] = bp_step_par[num_step][3];

    ioffset3 = (num_step - 1) * 3;
    cpxyz(mfoi, mst_org + ioffset3);
    cpxyz(hfoi, mst_orgH + ioffset3);

    ioffset9 = (num_step - 1) * 9;
    mst2orien(mst_orien, ioffset9, mfi);
    mst2orien(mst_orienH, ioffset9, hfi);
  }

  if (bz && bz_junction) {
    *str_type = -2; /* A-|B-/Z-DNA junction */
    if (z_step > 1)
      *str_type = -12;
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "Local base-pair step parameters\n");
  prt_stepstr(step_str, num_step, bphlx, 0, bp_step_par, fp);

  print_sep(fp, '*', 76);
  fprintf(fp, "Local base-pair helical parameters\n");
  prt_stepstr(step_str, num_step, bphlx, 1, bp_heli_par, fp);

  if (istep != 1) {
    print_sep(fp, '*', 76);
    fprintf(fp, "Step size [%ld]:  ", istep);
    fprintf(fp, "dinucleotide step classification based on Zp and ZpH will\n");
    fprintf(fp, "    be meaningless. Files <bp_helical.par>, <bp_step.par>, "
                "<hstacking.pdb>,\n");
    fprintf(fp, "    and <stacking.pdb> will also be corrupted.\n");
  }

  /* base-pair & step parameters for rebuilding */
  fstep = open_file(BPSTEP_FILE, "w");
  sprintf(str,
          "%4ld # base-pairs\n"
          "%4ld # ***local base-pair & step parameters***\n",
          num_bp, 0L);
  strcat(str, "#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening"
              "     Shift     Slide     Rise      Tilt      Roll      Twist\n");
  print_ds_rebuild_pars(bp_par, bp_step_par, num_bp, str, bp_seq, fstep);
  close_file(fstep);

  /* base-pair & helical parameters for rebuilding */
  fheli = open_file(HLXSTEP_FILE, "w");
  sprintf(str,
          "%4ld # base-pairs\n"
          "%4ld # ***local base-pair & helical parameters***\n",
          num_bp, 1L);
  strcat(str, "#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening"
              "     X-disp    Y-disp    h-Rise    Incl.     Tip     h-Twist\n");
  print_ds_rebuild_pars(bp_par, bp_heli_par, num_bp, str, bp_seq, fheli);
  close_file(fheli);

  /* for checking out */
  fchek = open_file(AUX_FILE, "w");
  fprintf(fchek, "Reference frame: Origins (Ox, Oy, Oz) followed by the"
                 " direction cosines of the\n"
                 "                 X- (Xx, Yx, Zx), Y- (Yx, Yx, Yx),"
                 " and Z- (Zx, Zx, Zx) axes\n");

  print_sep(fchek, '*', 89);
  fprintf(fchek, "Local base-pair reference frames\n");
  print_ref(bp_seq, num_bp, 1, bp_org, bp_orien, fchek);

  print_sep(fchek, '*', 89);
  fprintf(fchek, "Local middle reference frames\n");
  print_ref(bp_seq, num_bp - 1, 4, mst_org, mst_orien, fchek);

  print_sep(fchek, '*', 89);
  fprintf(fchek, "Local middle helical reference frames\n");
  print_ref(bp_seq, num_bp - 1, 4, mst_orgH, mst_orienH, fchek);

  print_sep(fchek, '*', 89);
  fprintf(fchek, "Local strand I base reference frames\n");
  print_ref(bp_seq, num_bp, 2, org[1], orien[1], fchek);

  print_sep(fchek, '*', 89);
  fprintf(fchek, "Local strand II base reference frames\n");
  print_ref(bp_seq, num_bp, 3, org[2], orien[2], fchek);

  print_sep(fchek, '*', 76);
  fprintf(fchek, "Local strand I base step parameters\n");
  parvec2mtx(step_par[1], nbpm1, bp_step_par);
  print_par(bp_seq, num_bp, 3, 0, bp_step_par, fchek);

  print_sep(fchek, '*', 76);
  fprintf(fchek, "Local strand I base helical parameters\n");
  parvec2mtx(heli_par[1], nbpm1, bp_heli_par);
  print_par(bp_seq, num_bp, 3, 1, bp_heli_par, fchek);

  print_sep(fchek, '*', 76);
  fprintf(fchek, "Local strand II base step parameters\n");
  parvec2mtx(step_par[2], nbpm1, bp_step_par);
  print_par(bp_seq, num_bp, 4, 0, bp_step_par, fchek);

  print_sep(fchek, '*', 76);
  fprintf(fchek, "Local strand II base helical parameters\n");
  parvec2mtx(heli_par[2], nbpm1, bp_heli_par);
  print_par(bp_seq, num_bp, 4, 1, bp_heli_par, fchek);

  close_file(fchek);

  /* reference frame for reseting the structure */
  print_analyze_ref_frames(2, num_bp, bp_seq, bp_org, bp_orien, pair_num,
                           nt_info);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dmatrix(hfi, 1, 3, 1, 3);
  free_dmatrix(bp_step_par, 1, nbpm1, 1, 6);
  free_dmatrix(bp_heli_par, 1, nbpm1, 1, 6);
  free_dmatrix(bp_par, 1, num_bp, 1, 6);
  free_dvector(bp_org, 1, num_bp * 3);
  free_dvector(bp_orien, 1, num_bp * 9);
  free_cmatrix(step_str, 1, num_bp, 0, 5);
}

void bz_check(double **r1, double *o1, double **r2, double *o2, long bz,
              long *bz_junction, long *z_step) {
  double x1[4], y1_osx[4], z1[4], x2[4], y2[4], z2[4], dorg[4];
  long i;

  if (!bz) /* without bz option */
    return;

  for (i = 1; i <= 3; i++) {
    x1[i] = r1[i][1];
    y1_osx[i] = r1[i][2];
    z1[i] = r1[i][3];
    x2[i] = r2[i][1];
    y2[i] = r2[i][2];
    z2[i] = r2[i][3];
    dorg[i] = o2[i] - o1[i];
  }

  if (dot(x1, x2) < 0.0 && dot(z1, z2) < 0.0 && dot(y1_osx, y2) > 0.0) {
    if (dot(dorg, z1) > 0.0) { /* bp2 in Z-form */
      for (i = 1; i <= 3; i++) {
        r2[i][1] = -r2[i][1];
        r2[i][3] = -r2[i][3];
      }
    } else { /* bp1 in Z-form */
      for (i = 1; i <= 3; i++) {
        r1[i][1] = -r1[i][1];
        r1[i][3] = -r1[i][3];
      }
    }
    (*bz_junction)++;
    return;
  }

  if (dot(x1, x2) > 0.0 && dot(z1, z2) > 0.0 && dot(y1_osx, y2) > 0.0 &&
      dot(dorg, z1) < 0.0 && dot(dorg, z2) < 0.0) { /* Z-step */
    for (i = 1; i <= 3; i++) {
      r1[i][1] = -r1[i][1]; /* x1 */
      r1[i][3] = -r1[i][3]; /* z1 */
      r2[i][1] = -r2[i][1]; /* x2 */
      r2[i][3] = -r2[i][3]; /* z2 */
    }
    (*z_step)++;
  }
}

static double chk_twist_for_larger_segment(double sum, long num, long nbpm1,
                                           long *idx) {
  long i;
  double ave = 0.0;

  if (num) {
    ave = sum / num;
    for (i = 1; i < nbpm1; i++)
      if (idx[i] && idx[i + 1])
        break;
    if (i == nbpm1 && i != 1)
      ave = 0.0;
  }

  return ave;
}

/* Get mean base-pair step twist angle: excluding breaks and non_WC steps */
void get_mtwist(long nbpm1, long *bphlx, long *WC_info, double **twist_rise,
                double *twist_p, double *twist_n) {
  long i, num_p = 0, num_n = 0, *idx_p, *idx_n;
  double twist, tp = 0.0, tn = 0.0;

  idx_p = lvector(1, nbpm1);
  idx_n = lvector(1, nbpm1);

  for (i = 1; i <= nbpm1; i++) {
    if (!bphlx[i] && WC_info[i] == 2 && WC_info[i + 1] == 2) {
      twist = twist_rise[i][1];
      if (twist >= 0) {
        idx_p[i] = 1;
        num_p++;
        tp += twist;
      } else {
        idx_n[i] = 1;
        num_n++;
        tn += twist;
      }
    }
  }

  *twist_p = chk_twist_for_larger_segment(tp, num_p, nbpm1, idx_p);
  *twist_n = chk_twist_for_larger_segment(tn, num_n, nbpm1, idx_n);

  free_lvector(idx_p, 1, DUMMY);
  free_lvector(idx_n, 1, DUMMY);
}

/* Calculate and print out 3DNA recommended local parameters */
void get_parameters(long ds, long num_bp, char **bp_seq, double **orien,
                    double **org, long *WC_info, FILE *fp, double **twist_rise,
                    double *mst_orien, double *mst_org, double *mst_orienH,
                    double *mst_orgH, long *bphlx, long istart, long istep,
                    long bz, long *str_type, long **pair_num, char **nt_info) {
  double hfoi[4], mfoi[4], o1[4], o2[4];
  double **heli_par, **step_par;
  double **hfi, **mfi, **r1, **r2;
  long i, j, m, ioffset3, ioffset9, nbpm1;

  nbpm1 = num_bp - 1;

  /* step and helical parameters for each strand */
  step_par = dmatrix(1, ds, 1, nbpm1 * 6);
  heli_par = dmatrix(1, ds, 1, nbpm1 * 6);
  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);
  hfi = dmatrix(1, 3, 1, 3);

  for (i = 1; i <= ds; i++) {
    for (j = 1; j <= nbpm1; j++) {
      refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
      m = (j - 1) * 6;
      bpstep_par(r1, o1, r2, o2, step_par[i] + m, mfi, mfoi);
      helical_par(r1, o1, r2, o2, heli_par[i] + m, hfi, hfoi);

      if (ds == 1) { /* populate mst_orien/mst_org, mst_orienH/mst_orgH */
        twist_rise[j][1] = step_par[i][m + 6];
        twist_rise[j][2] = step_par[i][m + 3];

        ioffset3 = (j - 1) * 3;
        cpxyz(mfoi, mst_org + ioffset3);
        cpxyz(hfoi, mst_orgH + ioffset3);

        ioffset9 = (j - 1) * 9;
        mst2orien(mst_orien, ioffset9, mfi);
        mst2orien(mst_orienH, ioffset9, hfi);
      }
    }
  }

  if (ds == 1)
    single_helix(num_bp, bp_seq, step_par, heli_par, orien, org, fp, pair_num,
                 nt_info);
  else
    double_helix(num_bp, bp_seq, step_par, heli_par, orien, org, WC_info, fp,
                 twist_rise, mst_orien, mst_org, mst_orienH, mst_orgH, bphlx,
                 istart, istep, bz, str_type, pair_num, nt_info);

  free_dmatrix(step_par, 1, ds, 1, nbpm1 * 6);
  free_dmatrix(heli_par, 1, ds, 1, nbpm1 * 6);
  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dmatrix(hfi, 1, 3, 1, 3);
}

/* Change vector-wise parameters to a num-by-6 matrix */
void parvec2mtx(double *parvec, long num, double **parmtx) {
  long i, ioffset, j, nc = 6;

  for (i = 1; i <= num; i++) {
    ioffset = (i - 1) * nc;
    for (j = 1; j <= nc; j++)
      parmtx[i][j] = parvec[ioffset + j];
  }
}

/* Print parameters for the rebuilding of a single helix */
void print_ss_rebuild_pars(double **pars, long num_bp, char *str, char **bp_seq,
                           FILE *fp) {
  char *fmt = " %9.3f";
  long i, j;

  fprintf(fp, "%s", str);

  /* 1st base: 6 zeros */
  fprintf(fp, "%c ", bp_seq[1][1]);
  for (i = 1; i <= 6; i++)
    fprintf(fp, fmt, 0.0);
  fprintf(fp, "\n");

  for (i = 2; i <= num_bp; i++) {
    fprintf(fp, "%c ", bp_seq[1][i]);
    for (j = 1; j <= 6; j++)
      fprintf(fp, fmt, pars[i - 1][j]);
    fprintf(fp, "\n");
  }
}

/* Print parameters for the rebuilding of a duplex */
void print_ds_rebuild_pars(double **bp_par, double **step_par, long num_bp,
                           char *str, char **bp_seq, FILE *fp) {
  char *fmt = " %9.3f";
  long i, j;

  fprintf(fp, "%s", str);

  /* 1st base-pair: 6 zeros for step parameters */
  fprintf(fp, "%c%c%c ", bp_seq[1][1], bp_seq[0][1], bp_seq[2][1]);
  for (i = 1; i <= 6; i++)
    fprintf(fp, fmt, bp_par[1][i]);
  for (i = 1; i <= 6; i++)
    fprintf(fp, fmt, 0.0);
  fprintf(fp, "\n");

  for (i = 2; i <= num_bp; i++) {
    fprintf(fp, "%c%c%c ", bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
    for (j = 1; j <= 6; j++)
      fprintf(fp, fmt, bp_par[i][j]);
    for (j = 1; j <= 6; j++)
      fprintf(fp, fmt, step_par[i - 1][j]);
    fprintf(fp, "\n");
  }
}

/* Print local base and base-pair reference frames */
void print_ref(char **bp_seq, long num_item, long ich, double *org,
               double *orien, FILE *fp) {
  long i, ioffset3, ioffset9, j;

  if (!lval_in_range(ich, 1, 4))
    fatal("wrong option [%ld] for printing reference frames\n", ich);

  fprintf(fp, "                Ox      Oy      Oz     Xx    Xy    Xz"
              "   Yx    Yy    Yz    Zx    Zy    Zz\n");

  for (i = 1; i <= num_item; i++) {
    if (ich == 1) /* base-pair */
      fprintf(fp, "%4ld %c%c%c   ", i, bp_seq[1][i], bp_seq[0][i],
              bp_seq[2][i]);
    else if (ich == 4) /* step */
      fprintf(fp, "%4ld %c%c/%c%c ", i, bp_seq[1][i], bp_seq[1][i + 1],
              bp_seq[2][i + 1], bp_seq[2][i]);
    else /* base I or II */
      fprintf(fp, "%4ld %c     ", i, bp_seq[ich - 1][i]);

    ioffset3 = (i - 1) * 3;
    for (j = 1; j <= 3; j++)
      fprintf(fp, "%8.2f", org[ioffset3 + j]);

    ioffset9 = (i - 1) * 9;
    for (j = 1; j <= 9; j++)
      fprintf(fp, "%6.2f", orien[ioffset9 + j]);

    fprintf(fp, "\n");
  }
}

/* Write multiple dinucleotide structures w.r.t. middle frames */
void write_mst(long ds, long num_bp, long **pair_num, char **bp_seq,
               double *mst_orien, double *mst_org, long **seidx,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               double **xyz, char **Miscs, long **htm_water,
               double **twist_rise, char *strfile) {
  double rise, **mst, **xyz_residue;
  long i, inum, ioffset3, j, jr, k, m;
  long tnum_res1, tnum_res2, inum_base;
  long ivec[3], ivec0[BUF512], ivect[BUF512];
  FILE *fp;

  if (num_bp == 1)
    return;

  fp = open_file(strfile, "w");

  xyz_residue = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
  mst = dmatrix(1, 3, 1, 3);

  for (i = 1; i <= num_bp - 1; i++) {
    rise = twist_rise[i][2];
    inum = 0;
    ioffset3 = (i - 1) * 3;
    orien2mst(mst_orien, (i - 1) * 9, mst);

    fprintf(fp, "%6s    %4ld\n", "MODEL ", i);
    if (!dval_in_range(rise, 2.0, 9.0)) /* rise out of normal range */
      fprintf(fp,
              "REMARK    NB [%.2f] -- the following step is unlikely in"
              " stacking geometry!\n",
              rise);

    if (ds == 1)
      fprintf(fp, "REMARK    Section #%4.4ld %c/%c\n", i, bp_seq[1][i],
              bp_seq[1][i + 1]);
    else
      fprintf(fp, "REMARK    Section #%4.4ld %c%c/%c%c\n", i, bp_seq[1][i],
              bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);

    fprintf(fp, "REMARK    %s\n", Gvars.X3DNA_VER);

    ivec[1] = pair_num[1][i]; /* lower level */
    if (ds == 1)
      inum_base = 1;
    else {
      inum_base = 2;
      ivec[2] = pair_num[2][i];
    }
    tnum_res1 = attached_residues(inum_base, ivec, ivec0, seidx, xyz, htm_water,
                                  &Gvars.misc_pars);
    fprintf(fp, "REMARK    LOWER: %4ld", tnum_res1);
    for (j = 1; j <= tnum_res1; j++) {
      ivect[j] = ivec0[j];
      fprintf(fp, "%4ld", j);
    }
    fprintf(fp, "\n");

    ivec[1] = pair_num[1][i + 1]; /* upper level */
    if (ds == 1)
      inum_base = 1;
    else {
      inum_base = 2;
      ivec[2] = pair_num[2][i + 1];
    }
    tnum_res2 = attached_residues(inum_base, ivec, ivec0, seidx, xyz, htm_water,
                                  &Gvars.misc_pars);
    fprintf(fp, "REMARK    UPPER: %4ld", tnum_res2);
    for (j = 1; j <= tnum_res2; j++) {
      k = j + tnum_res1;
      ivect[k] = ivec0[j];
      fprintf(fp, "%4ld", k);
    }
    fprintf(fp, "\n");

    for (j = 1; j <= tnum_res1 + tnum_res2; j++) {
      jr = ivect[j];
      for (k = seidx[jr][1]; k <= seidx[jr][2]; k++) {
        m = k - seidx[jr][1] + 1;
        cpxyz(xyz[k], xyz_residue[m]);
      }
      change_xyz(0, &mst_org[ioffset3], mst, seidx[jr][2] - seidx[jr][1] + 1,
                 xyz_residue);
      pdb_record(seidx[jr][1], seidx[jr][2], &inum, 1, AtomName, ResName,
                 ChainID, ResSeq, xyz_residue, Miscs, fp);
    }
    fprintf(fp, "ENDMDL\n");
  }

  free_dmatrix(xyz_residue, 1, NUM_RESIDUE_ATOMS, 1, 3);
  free_dmatrix(mst, 1, 3, 1, 3);

  close_file(fp);
}

void print_xyzP(long parallel, long nbpm1, char **bp_seq, long **phos,
                double *mst_orien, double *mst_org, double **xyz, FILE *fp,
                char *title_str, double **aveP, long p_offset) {
  char *bstr = "    --- ", *fmt = "%8.2f%8.2f%8.2f";
  char str[BUF512], temp_str[BUF512];
  double P_mst1[4], P_mst2[4], temp[4];
  long i, ioffset3, ioffset9, ip1, ip2, j;

  print_sep(fp, '*', 76);
  fprintf(fp, "%s", title_str);

  for (i = 1; i <= nbpm1; i++) {
    sprintf(str, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
            bp_seq[2][i + 1], bp_seq[2][i]);

    ioffset3 = (i - 1) * 3;
    ioffset9 = (i - 1) * 9;

    ip1 = phos[1 + p_offset][i + 1]; /* strand I */
    if (ip1) {
      ddxyz(mst_org + ioffset3, xyz[ip1], temp);
      for (j = 1; j <= 3; j++)
        P_mst1[j] = dot(temp, &mst_orien[ioffset9 + (j - 1) * 3]);
      sprintf(temp_str, fmt, P_mst1[1], P_mst1[2], P_mst1[3]);
      strcat(str, temp_str);
    } else
      for (j = 1; j <= 3; j++)
        strcat(str, bstr);

    ip2 = (parallel) ? phos[2 + p_offset][i + 1]
                     : phos[2 + p_offset][i]; /* strand II */
    if (ip2) {
      ddxyz(mst_org + ioffset3, xyz[ip2], temp);
      for (j = 1; j <= 3; j++)
        P_mst2[j] = dot(temp, &mst_orien[ioffset9 + (j - 1) * 3]);
      if (!parallel) {          /* anti-parallel */
        P_mst2[2] = -P_mst2[2]; /* reverse y */
        P_mst2[3] = -P_mst2[3]; /* reverse z */
      }
      sprintf(temp_str, fmt, P_mst2[1], P_mst2[2], P_mst2[3]);
      strcat(str, temp_str);
    } else
      for (j = 1; j <= 3; j++)
        strcat(str, bstr);
    if (ip1 && ip2)
      avexyz(P_mst1, P_mst2, aveP[i]);

    fprintf(fp, "%s\n", str);
  }
}

/* Following Stephen Harvey:
 *          1  |  Zp - ZpA     chi - chiA |
 *   ABI = --- | ---------- + ----------- |
 *          2  | ZpB - ZpA    chiB - chiA |
 *   where: ZpA = 2.2, ZpB = -0.4
 *          chiA = -157 (203); chiB = -108 (252)
 *   ref: Table 1 of the A-DNA motif paper, JMB2000
 */
static double get_ABI(long idx, double Zp, double **chi_angle) {
  double ZpA = 2.2, ZpB = -0.4, chiA = 203, chiB = 252;
  double x11, x12, x21, x22, xave, ABI = EMPTY_NUMBER;
  double ZpAB, chiAB, tZp, tchi;

  x11 = chi_angle[1][idx];
  x12 = chi_angle[1][idx + 1];
  x21 = chi_angle[2][idx];
  x22 = chi_angle[2][idx + 1];

  if (x11 > EMPTY_CRITERION && x12 > EMPTY_CRITERION && x21 > EMPTY_CRITERION &&
      x22 > EMPTY_CRITERION) {
    x11 = get_chi360(x11);
    x12 = get_chi360(x12);
    x21 = get_chi360(x21);
    x22 = get_chi360(x22);
    if (in_trans(x11) && in_trans(x12) && in_trans(x21) && in_trans(x22)) {
      xave = (x11 + x12 + x21 + x22) / 4.0;
      ZpAB = ZpB - ZpA;
      chiAB = chiB - chiA;
      tZp = (Zp - ZpA) / ZpAB;
      tchi = (xave - chiA) / chiAB;
      ABI = 0.5 * (tZp + tchi);
    }
  }

  return ABI;
}

/* Calculate and print xyz coordinates of P atoms w.r.t. middle frame
   and middle helical frame. A dinucleotide step is classified as A-,
   B- or TA-like */
void print_PP(long parallel, double **twist_rise, long num_bp, char **bp_seq,
              long **phos, double *mst_orien, double *mst_org,
              double *mst_orienH, double *mst_orgH, double **xyz, long *WC_info,
              long *bphlx, long abi, long **chi, FILE *fp) {
  char *bstr = "    --- ", *fmt = "%8.2f%8.2f%8.2f";
  char str[BUF512], **step_info;
  double **aveH, **aveS, **chi_angle = NULL, *ABIval = NULL;
  long i, j, ip, nbpm1, p_offset = 0, ds = 2;
  long *strABT, *idx;

  FILE *fchek;

  nbpm1 = num_bp - 1;

  aveS = dmatrix(1, nbpm1, 1, 3);
  aveH = dmatrix(1, nbpm1, 1, 3);

  fchek = open_file(AUX_FILE, "a");

  /* added August 3, 2004: for O1P & O2P as well */
  p_offset += 2;
  sprintf(
      str,
      "xyz coordinates of O1P atoms w.r.t. the middle frame of each dimer\n\n");
  strcat(str, "    step       xI      yI      zI     xII     yII     zII\n");
  print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orien, mst_org, xyz, fchek, str,
             aveS, p_offset);

  sprintf(str, "xyz coordinates of O1P atoms w.r.t. the middle helix frame"
               " of each dimer\n\n");
  strcat(str, "    step      xIH     yIH     zIH     xIIH    yIIH    zIIH\n");
  print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orienH, mst_orgH, xyz, fchek,
             str, aveH, p_offset);

  p_offset += 2;
  sprintf(
      str,
      "xyz coordinates of O2P atoms w.r.t. the middle frame of each dimer\n\n");
  strcat(str, "    step       xI      yI      zI     xII     yII     zII\n");
  print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orien, mst_org, xyz, fchek, str,
             aveS, p_offset);

  sprintf(str, "xyz coordinates of O2P atoms w.r.t. the middle helix frame"
               " of each dimer\n\n");
  strcat(str, "    step      xIH     yIH     zIH     xIIH    yIIH    zIIH\n");
  print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orienH, mst_orgH, xyz, fchek,
             str, aveH, p_offset);
  /* -------------------------------------------------------------------- */

  p_offset = 0; /* column index offset */
  sprintf(
      str,
      "xyz coordinates of P atoms w.r.t. the middle frame of each dimer\n\n");
  strcat(str, "    step       xI      yI      zI     xII     yII     zII\n");
  print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orien, mst_org, xyz, fchek, str,
             aveS, p_offset);

  sprintf(str, "xyz coordinates of P atoms w.r.t. the middle helix frame"
               " of each dimer\n\n");
  strcat(str, "    step      xIH     yIH     zIH     xIIH    yIIH    zIIH\n");
  print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orienH, mst_orgH, xyz, fchek,
             str, aveH, p_offset);

  close_file(fchek);

  if (abi) {
    chi_angle = dmatrix(1, ds, 1, num_bp);
    init_dmatrix(chi_angle, 1, ds, 1, num_bp, EMPTY_NUMBER);
    get_chi_torsions(ds, num_bp, chi, xyz, chi_angle);
    ABIval = dvector(1, num_bp);
    init_dvector(ABIval, 1, num_bp, EMPTY_NUMBER);
  }

  /* classification of each dinucleotide step as A-, B-, TA- or others */
  step_info = cmatrix(1, nbpm1, 0, BUF512);
  strABT = lvector(1, nbpm1);
  idx = lvector(1, num_bp);
  print_sep(fp, '*', 76);
  fprintf(fp, "Classification of each dinucleotide step in"
              " a right-handed nucleic acid\n"
              "structure: A-like; B-like; TA-like, or other cases.\n\n");
  if (abi)
    fprintf(
        fp,
        "For definition of the A-B index (ABI), see Waters et al. (2016).\n"
        "``Transitions of Double-Stranded DNA Between the A- and B-Forms.''\n"
        "J. Phys. Chem. B, 120(33), pp8449–8456.\n\n");

  fprintf(fp,
          "    step       Xp      Yp      Zp     XpH     YpH     ZpH    Form");
  if (abi)
    fprintf(fp, "   ABI");
  fprintf(fp, "\n");

  for (i = 1; i <= nbpm1; i++) {
    sprintf(step_info[i], "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
            bp_seq[2][i + 1], bp_seq[2][i]);
    ip = phos[1][i + 1] && ((parallel) ? phos[2][i + 1] : phos[2][i]);
    if (ip && !bphlx[i]) {
      sprintf(str, fmt, aveS[i][1], aveS[i][2], aveS[i][3]);
      strcat(step_info[i], str);
      sprintf(str, fmt, aveH[i][1], aveH[i][2], aveH[i][3]);
      strcat(step_info[i], str);
      if (WC_info[i] && WC_info[i + 1] &&                /* WC geometry */
          dval_in_range(twist_rise[i][1], 10.0, 60.0) && /* right-handed */
          dval_in_range(twist_rise[i][2], 2.5, 5.5) &&   /* Rise in range */
          dval_in_range(aveS[i][1], -5.0, -0.5) &&       /* Xp */
          dval_in_range(aveS[i][2], 7.5, 10.0) &&        /* Yp */
          dval_in_range(aveS[i][3], -2.0, 3.5) &&        /* Zp */
          dval_in_range(aveH[i][1], -11.5, 2.5) &&       /* XpH */
          dval_in_range(aveH[i][2], 1.5, 10.0) &&        /* YpH */
          dval_in_range(aveH[i][3], -3.0, 9.0)) {        /* ZpH */
        if (aveS[i][3] >= 1.5)                           /* A-form */
          strABT[i] = 1;
        else if (aveH[i][3] >= 4.0) /* TA-form */
          strABT[i] = 3;
        else if (aveS[i][3] <= 0.5 && aveH[i][1] < 0.5) /* B-form */
          strABT[i] = 2; /* aveS[i][3] < 0.5 for C-DNA #47 */
        if (abi)
          ABIval[i] = get_ABI(i, aveS[i][3], chi_angle);
      }
    } else
      for (j = 1; j <= 7; j++)
        strcat(step_info[i], bstr);
  }
  idx[1] = BUF512;
  idx[num_bp] = BUF512;
  for (i = 2; i <= nbpm1; i++)
    idx[i] = strABT[i] - strABT[i - 1];
  for (i = 1; i <= nbpm1; i++) {
    if (strABT[i] && (num_bp == 2 || !idx[i] || !idx[i + 1])) {
      if (strABT[i] == 1)
        fprintf(fp, "%s     A", step_info[i]);
      else if (strABT[i] == 2)
        fprintf(fp, "%s     B", step_info[i]);
      else
        fprintf(fp, "%s  *TA*", step_info[i]);
    } else
      fprintf(fp, "%s      ", step_info[i]);
    if (abi && ABIval[i] > EMPTY_CRITERION)
      fprintf(fp, " %6.2f", ABIval[i]);
    fprintf(fp, "\n");
  }

  free_dmatrix(aveS, 1, nbpm1, 1, 3);
  free_dmatrix(aveH, 1, nbpm1, 1, 3);
  free_cmatrix(step_info, 1, nbpm1, 0, BUF512);
  free_lvector(strABT, 1, nbpm1);
  free_lvector(idx, 1, num_bp);

  if (abi) {
    free_dmatrix(chi_angle, 1, ds, 1, num_bp);
    free_dvector(ABIval, 1, num_bp);
  }

  /* write P xyz coordinates to "auxiliary.par" */
  print_axyz(num_bp, bp_seq, phos, "P", xyz);
}

/* Classify a structure as right-handed or left handed forms */
void str_classify(double twist_p, double twist_n, long str_type, long parallel,
                  long num_bp, FILE *fp) {
  long mhlx = 0;

  print_sep(fp, '*', 76);
  fprintf(fp, "Structure classification: \n\n");

  /* added on 2014-dec-02 */
  if (!twist_p && !twist_n && !parallel && str_type == 1) {
    fprintf(fp, "This is a left-handed Z-form structure (as in Z-RNA, 1t4x)\n");
    return;
  }

  if (num_bp == 1) {
    fprintf(fp, "This structure contains only one base-pair\n");
    return;
  }
  if (parallel) {
    fprintf(fp, "This is a parallel duplex structure\n");
    return;
  }
  if (!twist_p && !twist_n)
    return;

  if (str_type >= 10) {
    mhlx = 1;
    fprintf(fp, "This structure contains more than one helical regions\n");
  }

  if (twist_p && twist_n) {
    if (str_type == -2)
      fprintf(fp, "This structure seems to contain a B-Z junction\n");
    if (str_type == -12)
      fprintf(fp, "This structure contains a B-Z junction...\n\n");
    if (str_type < 0) {
      fprintf(fp, "Note that 3DNA determines this junction based *purely* on "
                  "the geometry of\n"
                  "your input structure. The junction is between a "
                  "right-handed fragment with\n"
                  "a left-handed one, not necessarily just B-DNA and Z-DNA, "
                  "even though most\n"
                  "likely that would be the case.\n");
      return;
    }
  }

  str_type %= 10;
  if (str_type == 2) {
    fprintf(fp, "This nucleic acid structure is *unusual*\n");
    return;
  }

  if (!twist_p && twist_n) {
    if (str_type == 1)
      fprintf(fp, "This is a left-handed Z-form structure\n");
    else if (!mhlx)
      fprintf(fp, "This is a left-handed W-form structure\n");
  }

  if (twist_p && !twist_n) {
    fprintf(fp, "This is a right-handed ");
    if (str_type == 1)
      fprintf(fp, "unknown R-form structure\n");
    else
      fprintf(fp, "nucleic acid structure\n");
  }
}

/* Calculate helix radius */
double a_hlxdist(long idx, double **xyz, double *hlx_axis, double *hlx_pos) {
  double temp, d[4];
  long i;

  if (idx) {
    ddxyz(hlx_pos, xyz[idx], d);
    temp = dot(d, hlx_axis);
    for (i = 1; i <= 3; i++)
      d[i] -= temp * hlx_axis[i];
    return veclen(d);
  } else
    return EMPTY_NUMBER;
}

/* Print the radius from P, O4' & C1' atoms to the local helical axis */
void print_radius(char **bp_seq, long nbpm1, long ich, double **p_radius,
                  double **o4_radius, double **c1_radius, long *bphlx,
                  FILE *fp) {
  char *bstr = "      ----", *fmt = "%10.2f";
  char str[BUF512];
  long i, j, ik;

  if (!lval_in_range(ich, 1, 3))
    fatal("wrong option [%ld] for printing helix radius\n", ich);

  if (ich == 1) { /* duplex */
    fprintf(fp, "                        Strand I"
                "                      Strand II\n"
                "     step         P        O4'       C1'"
                "        P        O4'        C1'\n");
    for (i = 1; i <= nbpm1; i++) {
      sprintf(str, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
              bp_seq[2][i + 1], bp_seq[2][i]);
      if (bphlx[i])
        for (j = 1; j <= 6; j++)
          strcat(str, bstr);
      else {
        for (j = 1; j <= 2; j++) {
          parcat(str, p_radius[j][i], fmt, bstr);
          parcat(str, o4_radius[j][i], fmt, bstr);
          parcat(str, c1_radius[j][i], fmt, bstr);
        }
      }
      fprintf(fp, "%s\n", str);
    }

  } else { /* single strand */
    fprintf(fp, "    step          P        O4'       C1'\n");
    ik = ich - 1;
    for (i = 1; i <= nbpm1; i++) {
      sprintf(str, "%4ld  %c/%c ", i, bp_seq[ik][i], bp_seq[ik][i + 1]);
      parcat(str, p_radius[ik][i], fmt, bstr);
      parcat(str, o4_radius[ik][i], fmt, bstr);
      parcat(str, c1_radius[ik][i], fmt, bstr);
      fprintf(fp, "%s\n", str);
    }
  }
}

/* Get radius from P, O4' and C1' to the helical axis */
void helix_radius(long ds, long num_bp, char **bp_seq, double **orien,
                  double **org, long **phos, long **chi, double **xyz,
                  long *bphlx, FILE *fp) {
  double temp1, temp2;
  double hx[4], morg[4], o1[44], o2[4], pars[7];
  double **mst, **r1, **r2, **c1_radius, **o4_radius, **p_radius;
  double *bp_orien, *bp_org, **c1BP_radius, **o4BP_radius, **pBP_radius;
  long i, ioffset, j, k, nbpm1, pn;
  FILE *fchek;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mst = dmatrix(1, 3, 1, 3);

  p_radius = dmatrix(1, ds, 1, nbpm1);
  o4_radius = dmatrix(1, ds, 1, nbpm1); /* two O4' atom per step */
  c1_radius = dmatrix(1, ds, 1, nbpm1); /* two C1' atom per step */

  for (i = 1; i <= ds; i++)
    for (j = 1; j <= nbpm1; j++) {
      refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
      helical_par(r1, o1, r2, o2, pars, mst, morg);
      for (k = 1; k <= 3; k++)
        hx[k] = mst[k][3];

      pn = (i == 1) ? j + 1 : j; /* P index: +1 for I */
      p_radius[i][j] = a_hlxdist(phos[i][pn], xyz, hx, morg);

      ioffset = (j - 1) * 4;

      temp1 = a_hlxdist(chi[i][ioffset + 1], xyz, hx, morg);
      temp2 = a_hlxdist(chi[i][ioffset + 5], xyz, hx, morg);
      o4_radius[i][j] = 0.5 * (temp1 + temp2);
      temp1 = a_hlxdist(chi[i][ioffset + 2], xyz, hx, morg);
      temp2 = a_hlxdist(chi[i][ioffset + 6], xyz, hx, morg);
      c1_radius[i][j] = 0.5 * (temp1 + temp2);
    }

  print_sep(fp, '*', 76);
  fprintf(fp, "Helix radius (radial displacement of P, O4', and C1'"
              " atoms in local helix\n   frame of each dimer)\n\n");

  if (ds == 1)
    print_radius(bp_seq, nbpm1, 2, p_radius, o4_radius, c1_radius, bphlx, fp);
  else {
    bp_org = dvector(1, num_bp * 3);
    bp_orien = dvector(1, num_bp * 9);

    pBP_radius = dmatrix(1, ds, 1, nbpm1);
    o4BP_radius = dmatrix(1, ds, 1, nbpm1);
    c1BP_radius = dmatrix(1, ds, 1, nbpm1);

    for (i = 1; i <= num_bp; i++) {
      refs_right_left(i, orien, org, r1, o1, r2, o2);
      bpstep_par(r1, o1, r2, o2, pars, mst, morg);

      cpxyz(morg, bp_org + (i - 1) * 3);
      mst2orien(bp_orien, (i - 1) * 9, mst);
    }

    for (i = 1; i <= nbpm1; i++) {
      refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
      helical_par(r1, o1, r2, o2, pars, mst, morg);
      for (k = 1; k <= 3; k++)
        hx[k] = mst[k][3];

      ioffset = (i - 1) * 4;

      for (j = 1; j <= ds; j++) {
        pn = (j == 1) ? i + 1 : i; /* P index: +1 for I */

        pBP_radius[j][i] = a_hlxdist(phos[j][pn], xyz, hx, morg);

        temp1 = a_hlxdist(chi[j][ioffset + 1], xyz, hx, morg);
        temp2 = a_hlxdist(chi[j][ioffset + 5], xyz, hx, morg);
        o4BP_radius[j][i] = 0.5 * (temp1 + temp2);

        temp1 = a_hlxdist(chi[j][ioffset + 2], xyz, hx, morg);
        temp2 = a_hlxdist(chi[j][ioffset + 6], xyz, hx, morg);
        c1BP_radius[j][i] = 0.5 * (temp1 + temp2);
      }
    }
    print_radius(bp_seq, nbpm1, 1, pBP_radius, o4BP_radius, c1BP_radius, bphlx,
                 fp);

    fchek = open_file(AUX_FILE, "a");

    print_sep(fchek, '*', 76);
    fprintf(fchek, "Strand I helix radius\n\n");
    print_radius(bp_seq, nbpm1, 2, p_radius, o4_radius, c1_radius, bphlx,
                 fchek);

    print_sep(fchek, '*', 76);
    fprintf(fchek, "Strand II helix radius\n\n");
    print_radius(bp_seq, nbpm1, 3, p_radius, o4_radius, c1_radius, bphlx,
                 fchek);

    close_file(fchek);

    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_dmatrix(pBP_radius, 1, ds, 1, nbpm1);
    free_dmatrix(o4BP_radius, 1, ds, 1, nbpm1);
    free_dmatrix(c1BP_radius, 1, ds, 1, nbpm1);
  }

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mst, 1, 3, 1, 3);
  free_dmatrix(p_radius, 1, ds, 1, nbpm1);
  free_dmatrix(o4_radius, 1, ds, 1, nbpm1);
  free_dmatrix(c1_radius, 1, ds, 1, nbpm1);
}

/* Print local base/base-pair helical reference frames */
void print_shlx(char **bp_seq, long nbpm1, long ich, double *shlx_orien,
                double *shlx_org, FILE *fp) {
  long i, ioffset6, ioffset18, j, k;

  if (!lval_in_range(ich, 1, 3))
    fatal("wrong option [%ld] for printing helix axis\n", ich);

  fprintf(fp, "               Ox      Oy      Oz     Xx    Xy    Xz"
              "   Yx    Yy    Yz    Zx    Zy    Zz\n");

  k = ich - 1;
  for (i = 1; i <= nbpm1; i++) {
    if (ich == 1)
      fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
              bp_seq[2][i + 1], bp_seq[2][i]);
    else
      fprintf(fp, "%4ld  %c/%c ", i, bp_seq[k][i], bp_seq[k][i + 1]);

    ioffset6 = (i - 1) * 6;
    ioffset18 = (i - 1) * 18;
    for (j = 1; j <= 3; j++)
      fprintf(fp, "%8.2f", shlx_org[ioffset6 + j]);
    for (j = 1; j <= 9; j++)
      fprintf(fp, "%6.2f", shlx_orien[ioffset18 + j]);
    fprintf(fp, "\n          ");

    ioffset6 += 3;
    ioffset18 += 9;
    for (j = 1; j <= 3; j++)
      fprintf(fp, "%8.2f", shlx_org[ioffset6 + j]);
    for (j = 1; j <= 9; j++)
      fprintf(fp, "%6.2f", shlx_orien[ioffset18 + j]);
    fprintf(fp, "\n");
  }
}

/* Get the local helical axis for bending analysis */
void get_helix_axis(long ds, long num_bp, char **bp_seq, double **orien,
                    double **org, long *bphlx, FILE *fp) {
  char *bstr = "      ----", *fmt = "%10.2f";
  double hf_twist, hf_rise;
  double hx[4], morg[4], o1[4], o2[4], pars[7];
  double **hlx_org, **hlx_orien, **mst, **r1, **r2, **temp1, **temp2;
  double *bp_hlx_org, *bp_hlx_orien, *bp_org, *bp_orien;
  long i, ioffset, j, k, nbpm1;
  FILE *fchek;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mst = dmatrix(1, 3, 1, 3);
  temp1 = dmatrix(1, 3, 1, 3);
  temp2 = dmatrix(1, 3, 1, 3);

  hlx_orien = dmatrix(1, ds, 1, nbpm1 * 18);
  hlx_org = dmatrix(1, ds, 1, nbpm1 * 6);

  for (i = 1; i <= ds; i++)
    for (j = 1; j <= nbpm1; j++) {
      refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
      helical_par(r1, o1, r2, o2, pars, mst, morg);
      for (k = 1; k <= 3; k++)
        hx[k] = mst[k][3];

      hf_twist = 0.5 * pars[6];
      hf_rise = 0.5 * pars[3];

      ioffset = (j - 1) * 18;
      arb_rotation(hx, -hf_twist, temp1);
      multi_matrix(temp1, 3, 3, mst, 3, 3, temp2);
      mst2orien(hlx_orien[i], ioffset, temp2);

      ioffset += 9;
      arb_rotation(hx, hf_twist, temp1);
      multi_matrix(temp1, 3, 3, mst, 3, 3, temp2);
      mst2orien(hlx_orien[i], ioffset, temp2);

      ioffset = (j - 1) * 6;
      for (k = 1; k <= 3; k++) {
        hlx_org[i][ioffset + k] = morg[k] - hf_rise * hx[k];
        hlx_org[i][ioffset + k + 3] = morg[k] + hf_rise * hx[k];
      }
    }

  fchek = open_file(AUX_FILE, "a");

  if (ds == 1) {
    print_sep(fp, '*', 76);
    fprintf(fp, "Position (Px, Py, Pz) and local helical axis vector"
                " (Hx, Hy, Hz)\n\n");
    fprintf(fp, "     step       Px        Py        Pz"
                "        Hx        Hy        Hz\n");

    for (i = 1; i <= nbpm1; i++) {
      ioffset = (i - 1) * 6;
      avexyz(hlx_org[ds] + ioffset, hlx_org[ds] + ioffset + 3, morg);
      cpxyz(hlx_orien[ds] + (i - 1) * 18 + 6, hx);

      fprintf(fp, "%4ld  %c/%c ", i, bp_seq[ds][i], bp_seq[ds][i + 1]);
      for (j = 1; j <= 3; j++)
        fprintf(fp, fmt, morg[j]);
      for (j = 1; j <= 3; j++)
        fprintf(fp, fmt, hx[j]);
      fprintf(fp, "\n");
    }

    print_sep(fchek, '*', 76);
    fprintf(fchek, "Helix axis\n\n");
    print_shlx(bp_seq, nbpm1, 2, hlx_orien[ds], hlx_org[ds], fchek);

  } else { /* duplex */
    bp_orien = dvector(1, num_bp * 9);
    bp_org = dvector(1, num_bp * 3);
    bp_hlx_orien = dvector(1, nbpm1 * 18);
    bp_hlx_org = dvector(1, nbpm1 * 6);

    for (i = 1; i <= num_bp; i++) {
      refs_right_left(i, orien, org, r1, o1, r2, o2);
      bpstep_par(r1, o1, r2, o2, pars, mst, morg);

      cpxyz(morg, bp_org + (i - 1) * 3);
      mst2orien(bp_orien, (i - 1) * 9, mst);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Position (Px, Py, Pz) and local helical axis vector"
                " (Hx, Hy, Hz)\n         for each dinucleotide step\n\n");
    fprintf(fp, "     step       Px        Py        Pz"
                "        Hx        Hy        Hz\n");

    for (i = 1; i <= nbpm1; i++) {
      refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
      helical_par(r1, o1, r2, o2, pars, mst, morg);
      for (j = 1; j <= 3; j++)
        hx[j] = mst[j][3];

      hf_twist = 0.5 * pars[6];
      hf_rise = 0.5 * pars[3];

      ioffset = (i - 1) * 18;
      arb_rotation(hx, -hf_twist, temp1);
      multi_matrix(temp1, 3, 3, mst, 3, 3, temp2);
      mst2orien(bp_hlx_orien, ioffset, temp2);

      ioffset += 9;
      arb_rotation(hx, hf_twist, temp1);
      multi_matrix(temp1, 3, 3, mst, 3, 3, temp2);
      mst2orien(bp_hlx_orien, ioffset, temp2);

      ioffset = (i - 1) * 6;
      for (j = 1; j <= 3; j++) {
        bp_hlx_org[ioffset + j] = morg[j] - hf_rise * hx[j];
        bp_hlx_org[ioffset + j + 3] = morg[j] + hf_rise * hx[j];
      }

      fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
              bp_seq[2][i + 1], bp_seq[2][i]);
      if (bphlx[i])
        for (j = 1; j <= 6; j++)
          fprintf(fp, "%s", bstr);
      else {
        for (j = 1; j <= 3; j++)
          fprintf(fp, fmt, morg[j]);
        for (j = 1; j <= 3; j++)
          fprintf(fp, fmt, hx[j]);
      }
      fprintf(fp, "\n");
    }

    print_sep(fchek, '*', 88);
    fprintf(fchek, "Base-pair helix axis\n\n");
    print_shlx(bp_seq, nbpm1, 1, bp_hlx_orien, bp_hlx_org, fchek);

    print_sep(fchek, '*', 88);
    fprintf(fchek, "Strand I helix axis\n\n");
    print_shlx(bp_seq, nbpm1, 2, hlx_orien[1], hlx_org[1], fchek);

    print_sep(fchek, '*', 88);
    fprintf(fchek, "Strand II helix axis\n\n");
    print_shlx(bp_seq, nbpm1, 3, hlx_orien[2], hlx_org[2], fchek);

    free_dvector(bp_orien, 1, num_bp * 9);
    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_hlx_orien, 1, nbpm1 * 18);
    free_dvector(bp_hlx_org, 1, nbpm1 * 6);
  }

  close_file(fchek);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mst, 1, 3, 1, 3);
  free_dmatrix(temp1, 1, 3, 1, 3);
  free_dmatrix(temp2, 1, 3, 1, 3);
  free_dmatrix(hlx_orien, 1, ds, 1, nbpm1 * 18);
  free_dmatrix(hlx_org, 1, ds, 1, nbpm1 * 6);
}

/* Get the helical axis based on C1* and RN9/YN1 atoms */
void get_axis(long nvec, long **idx, long num, double **xyz, long nb, long *C1b,
              long *C1e, double *std_rise, double *hrise, double *haxis,
              double *hstart, double *hend) {
  double tb, te, hinge[4], org_xyz[4], t2[3];
  double *g, *drise;
  double **dxy, **dxyT, **rotmat, **vxyz, **xyzH;
  double **dd, **inv_dd;
  long i, j;

  *hrise = EMPTY_NUMBER;
  if (nvec < 3) /* too few vectors to define a helix */
    return;

  /* find helical axis and rise */
  vxyz = dmatrix(1, nvec, 1, 3);
  drise = dvector(1, nvec);
  for (i = 1; i <= nvec; i++)
    ddxyz(xyz[idx[i][1]], xyz[idx[i][2]], vxyz[i]);
  ls_plane(vxyz, nvec, haxis, hinge, hrise, drise);
  if (*hrise < 0.0) {
    *hrise = -*hrise;
    negate_xyz(haxis);
  }
  *std_rise = std_dvector(drise, nvec);

  /* align haxis to global z-axis */
  rotmat = dmatrix(1, 3, 1, 3);
  xyzH = dmatrix(1, num, 1, 3);
  align2zaxis(num, haxis, rotmat, xyz, xyzH);

  /* locate xy-coordinate the helix passes through */
  dxy = dmatrix(1, nvec, 1, 2);
  dxyT = dmatrix(1, 2, 1, nvec);
  g = dvector(1, nvec);
  for (i = 1; i <= nvec; i++) {
    g[i] = 0.0;
    for (j = 1; j <= 2; j++) {
      tb = xyzH[idx[i][1]][j];
      te = xyzH[idx[i][2]][j];
      dxy[i][j] = 2.0 * (te - tb);
      g[i] += te * te - tb * tb;
    }
  }

  dd = dmatrix(1, 2, 1, 2);
  inv_dd = dmatrix(1, 2, 1, 2);
  multi_vec_matrix(g, nvec, dxy, nvec, 2, t2);
  transpose_matrix(dxy, nvec, 2, dxyT);
  multi_matrix(dxyT, 2, nvec, dxy, nvec, 2, dd);
  dinverse(dd, 2, inv_dd);
  multi_vec_matrix(t2, 2, inv_dd, 2, 2, org_xyz);

  tb = 0.0; /* ave. z-coordinate of 1st base set C1* atoms */
  j = 0;
  for (i = 1; i <= nb; i++) {
    if (C1b[i]) {
      j++;
      tb += xyzH[C1b[i]][3];
    }
  }
  org_xyz[3] = tb / j;
  multi_vec_matrix(org_xyz, 3, rotmat, 3, 3, hstart);

  te = 0.0; /* ave. z-coordinate of last base set C1* atoms */
  j = 0;
  for (i = 1; i <= nb; i++) {
    if (C1e[i]) {
      j++;
      te += xyzH[C1e[i]][3];
    }
  }
  org_xyz[3] = te / j;
  multi_vec_matrix(org_xyz, 3, rotmat, 3, 3, hend);

  free_dmatrix(vxyz, 1, nvec, 1, 3);
  free_dvector(drise, 1, nvec);
  free_dmatrix(rotmat, 1, 3, 1, 3);
  free_dmatrix(xyzH, 1, num, 1, 3);
  free_dmatrix(dxy, 1, nvec, 1, 2);
  free_dmatrix(dxyT, 1, 2, 1, nvec);
  free_dvector(g, 1, nvec);
  free_dmatrix(dd, 1, 2, 1, 2);
  free_dmatrix(inv_dd, 1, 2, 1, 2);
}

/* Print out P, O4' and C1' atom radius defined cylinder in Raster3D format */
void print_poc_r3d(double *rave, double *hstart, double *hend) {
  char label_style[BUF512];
  double width3[4], hb_col[5], **atom_col, **base_col;
  long i;
  FILE *fpr3d;

  atom_col = dmatrix(0, NATOMCOL, 1, 3);
  base_col = dmatrix(0, NBASECOL, 1, 3);
  get_r3dpars(base_col, hb_col, width3, atom_col, label_style);

  fpr3d = open_file(POC_FILE, "w");
  fprintf(fpr3d,
          "###\n### Linear global helical axis if not strongly curved\n");

  rave[3] = width3[1];
  for (i = 0; i < 4;
       i++) /* 0 for P; 1 for O, 2 for C, and 3 for thinner line */
    r3d_rod((i != 3) ? 99L : 5L, hstart, hend, rave[i], hb_col, fpr3d);
  close_file(fpr3d);

  free_dmatrix(atom_col, 0, NATOMCOL, 1, 3);
  free_dmatrix(base_col, 0, NBASECOL, 1, 3);
}

/* Get the radius of P, O4' and C1' atoms from the global helical axis */
static void get_poc_radius(long ds, long num_bp, long **phos, long **chi,
                           double **xyz, double *haxis, double *hstart,
                           double *hend, double *rave, double *rstd) {
  long i, ik, ioffset, j, k, num_ple, poc_num[3];
  double **poc_radius;

  for (i = 0; i <= 2; i++) {
    poc_num[i] = 0;
    rave[i] = rstd[i] = 0.0;
  }

  num_ple = 2 * num_bp; /* maximum possible number */
  poc_radius = dmatrix(0, 2, 1, num_ple);
  for (i = 1; i <= ds; i++)
    for (j = 1; j <= num_bp; j++) {
      k = 0; /* P */
      ik = phos[i][j];
      if (ik)
        poc_radius[k][++poc_num[k]] = a_hlxdist(ik, xyz, haxis, hstart);
      ioffset = (j - 1) * 4;
      k = 1; /* O4 */
      ik = chi[i][ioffset + 1];
      if (ik)
        poc_radius[k][++poc_num[k]] = a_hlxdist(ik, xyz, haxis, hstart);
      k = 2; /* C1 */
      ik = chi[i][ioffset + 2];
      if (ik)
        poc_radius[k][++poc_num[k]] = a_hlxdist(ik, xyz, haxis, hstart);
    }

  for (i = 0; i <= 2; i++)
    if (poc_num[i]) {
      rave[i] = ave_dvector(poc_radius[i], poc_num[i]);
      rstd[i] = std_dvector(poc_radius[i], poc_num[i]);
    }
  print_poc_r3d(rave, hstart, hend);

  free_dmatrix(poc_radius, 0, 2, 1, num_ple);
}

/* Global helical parameters based on C1'-C1' vector for a duplex */
static void global_c1_c1_par(long num_bp, char **bp_seq, long **chi,
                             double **xyz, double *haxis, double *hstart,
                             FILE *fp) {
  char *bstr = "      --- ", *fmt = "%10.2f", **c1_hpar;
  double dtmp, dd[4], org_xyz[4], **dc1, **mc1;
  long i, ia, ib, ioffset, j;

  fprintf(fp, "\nGlobal parameters based on C1'-C1' vectors:\n\n");
  fprintf(
      fp,
      "disp.: displacement of the middle C1'-C1' point from the helix\n"
      "angle: inclination between C1'-C1' vector and helix (subtracted from "
      "90)\n"
      "twist: helical twist angle between consecutive C1'-C1' vectors\n"
      "rise:  helical rise by projection of the vector connecting consecutive\n"
      "       C1'-C1' middle points onto the helical axis\n\n");

  c1_hpar = cmatrix(1, num_bp, 0, 50);
  dc1 = dmatrix(1, num_bp, 1, 3);
  mc1 = dmatrix(1, num_bp, 1, 3);
  for (i = 1; i <= num_bp; i++) { /* displacement and angle */
    sprintf(c1_hpar[i], "%4ld %c%c%c", i, bp_seq[1][i], bp_seq[0][i],
            bp_seq[2][i]);
    ioffset = (i - 1) * 4;
    ia = chi[1][ioffset + 2];
    ib = chi[2][ioffset + 2];
    if (ia && ib) {
      ddxyz(xyz[ib], xyz[ia], dc1[i]);
      avexyz(xyz[ia], xyz[ib], mc1[i]);
      ddxyz(hstart, mc1[i], dd);
      dtmp = dot(dd, haxis);
      for (j = 1; j <= 3; j++)
        org_xyz[j] = dd[j] - dtmp * haxis[j];
      parcat(c1_hpar[i], veclen(org_xyz), fmt, bstr);
      dtmp = 90.0 - magang(dc1[i], haxis);
      parcat(c1_hpar[i], dtmp, fmt, bstr);
    } else {
      parcat(c1_hpar[i], EMPTY_NUMBER, fmt, bstr);
      parcat(c1_hpar[i], EMPTY_NUMBER, fmt, bstr);
    }
  }

  for (i = 1; i <= num_bp - 1; i++) { /* twist and rise */
    if (strstr(c1_hpar[i], bstr) == NULL &&
        strstr(c1_hpar[i + 1], bstr) == NULL) {
      dtmp = vec_ang(dc1[i], dc1[i + 1], haxis);
      parcat(c1_hpar[i], dtmp, fmt, bstr);
      ddxyz(mc1[i], mc1[i + 1], dd);
      parcat(c1_hpar[i], dot(dd, haxis), fmt, bstr);
    } else {
      parcat(c1_hpar[i], EMPTY_NUMBER, fmt, bstr);
      parcat(c1_hpar[i], EMPTY_NUMBER, fmt, bstr);
    }
  }
  parcat(c1_hpar[num_bp], EMPTY_NUMBER, fmt, bstr);
  parcat(c1_hpar[num_bp], EMPTY_NUMBER, fmt, bstr);

  fprintf(fp, "     bp       disp.    angle     twist      rise\n");
  for (i = 1; i <= num_bp; i++)
    fprintf(fp, "%s\n", c1_hpar[i]);

  free_cmatrix(c1_hpar, 1, num_bp, 0, 50);
  free_dmatrix(dc1, 1, num_bp, 1, 3);
  free_dmatrix(mc1, 1, num_bp, 1, 3);
}

/* Structural analysis from a global prospective of the backbone:
   measuring curvature, and bending angle */
void global_analysis(long ds, long num_bp, long num, char **bp_seq, long **chi,
                     long **phos, double **xyz, FILE *fp) {
  double hrise, std_rise, haxis[4], hstart[4], hend[4], rave[4], rstd[3];
  long nvec, C1b[MBASES], C1e[MBASES], **idx;

  idx = lmatrix(1, 4 * num_bp, 1, 2); /* beginning & end index */
  get_CNidx(ds, num_bp, chi, idx, &nvec, C1b, C1e);
  get_axis(nvec, idx, num, xyz, ds, C1b, C1e, &std_rise, &hrise, haxis, hstart,
           hend);
  free_lmatrix(idx, 1, 4 * num_bp, 1, 2);

  if (hrise < EMPTY_CRITERION)
    return; /* no helix defined */

  print_sep(fp, '*', 76);
  fprintf(fp, "Global linear helical axis defined by equivalent C1'"
              " and RN9/YN1 atom pairs\n");
  fprintf(fp, "Deviation from regular linear helix: %.2f(%.2f)\n", hrise,
          std_rise);

  if (std_rise > Gvars.misc_pars.std_curved)
    return;

  get_poc_radius(ds, num_bp, phos, chi, xyz, haxis, hstart, hend, rave, rstd);

  fprintf(fp, "Helix:  %9.4f %9.4f %9.4f\n", haxis[1], haxis[2], haxis[3]);
  fprintf(fp, "HETATM 9998  XS    X X 999    %8.3f%8.3f%8.3f\n", hstart[1],
          hstart[2], hstart[3]);
  fprintf(fp, "HETATM 9999  XE    X X 999    %8.3f%8.3f%8.3f\n", hend[1],
          hend[2], hend[3]);
  fprintf(fp, "Average and standard deviation of helix radius:\n");
  fprintf(fp, "      P: %.2f(%.2f), O4': %.2f(%.2f),  C1': %.2f(%.2f)\n",
          rave[0], rstd[0], rave[1], rstd[1], rave[2], rstd[2]);

  if (ds == 1)
    return;

  global_c1_c1_par(num_bp, bp_seq, chi, xyz, haxis, hstart, fp);
}

/* Check if two nucleotides along single strand are reasonably stacked
 * [oct-12-2007] */
static long with_ss_overlap(long istep, long r1, long r2, double **orien,
                            long **ring_atom, double **xyz) {
  long i, ioffset, j, n1, n2;
  double sdist = XBIG, oave[4] = {EMPTY_NUMBER, 0.0, 0.0, 0.0};
  double **oxyz1, **oxyz2;

  ioffset = (istep - 1) * 9;

  if (z1_z2_angle_in_0_to_90(&orien[1][ioffset + 6], &orien[1][ioffset + 15]) >
      65.0)
    return 0;

  oxyz1 = dmatrix(1, 9, 1, 3);
  oxyz2 = dmatrix(1, 9, 1, 3);

  n1 = ratom_xyz(ring_atom[r1], 0, xyz, oave,
                 oxyz1); /* including exocyclic atoms */
  n2 = ratom_xyz(ring_atom[r2], 0, xyz, oave, oxyz2);

  for (i = 1; i <= n1; i++)
    for (j = 1; j <= n2; j++)
      sdist = dval_min(sdist, p1p2_dist(oxyz1[i], oxyz2[j]));

  free_dmatrix(oxyz1, 1, 9, 1, 3);
  free_dmatrix(oxyz2, 1, 9, 1, 3);

  if (sdist > 4.5)
    return 0;

  return 1;
}

static long with_stacking(long ds, long istep, long **pair_num, double **orien,
                          long **ring_atom, double **xyz) {
  long i1, i2, k1, k2;

  i1 = pair_num[1][istep];
  i2 = pair_num[1][istep + 1];

  if (ds == 1)
    return with_ss_overlap(istep, i1, i2, orien, ring_atom, xyz);

  k1 = pair_num[ds][istep];
  k2 = pair_num[ds][istep + 1];

  return with_ss_overlap(istep, i1, i2, orien, ring_atom, xyz) ||
         with_ss_overlap(istep, k2, k1, orien, ring_atom, xyz) ||
         with_ss_overlap(istep, i1, k2, orien, ring_atom, xyz) ||
         with_ss_overlap(istep, i2, k1, orien, ring_atom, xyz);
}

/* Get overlap area between neighbor bases (ring + first order of shell) */
void base_overlap(long ds, long num_bp, long num, long num_residue,
                  long **pair_num, long *bRY, char **bp_seq, long **seidx,
                  char **AtomName, double **xyz, long *idx, double **orien,
                  double **org, FILE *fp) {
  char *fmt = " %5.2f(%5.2f)";
  long i, j, k, r1, r2;
  long **ring_atom;
  double oave[4], zave[4], *olarea;

  /* 1-9 ring atom index, 10 # of ring atoms, 11-19 first level */
  ring_atom = lmatrix(1, num_residue, 1, 19);
  ring_oidx(num, num_residue, bRY, seidx, AtomName, xyz, idx, ring_atom);

  /* i1-i2(1), i1-j2(2), j1-i2(3), j1-j2(4) ==> sum(5): ring + 1st level */
  olarea = dvector(1, 10);

  print_sep(fp, '*', 76);
  fprintf(
      fp,
      "Overlap area in Angstrom^2 between polygons defined by atoms on"
      " successive\nbases. Polygons projected in the mean plane of the designed"
      " base-pair step.\n\n");
  fprintf(
      fp,
      "Values in parentheses measure the overlap of base ring atoms only."
      " Those\noutside parentheses include exocyclic atoms on the ring. Intra-"
      " and\ninter-strand overlap is designated according to the following"
      " diagram:\n\n");
  if (ds == 2) {
    fprintf(fp, "                    i2  3'      5' j2\n"
                "                       /|\\      |\n"
                "                        |       |\n"
                "               Strand I |       | II\n"
                "                        |       |\n"
                "                        |       |\n"
                "                        |      \\|/\n"
                "                    i1  5'      3' j1\n\n");
    fprintf(fp, "     step      i1-i2        i1-j2        j1-i2        j1-j2"
                "        sum\n");
  } else {
    fprintf(fp, "                    i2  3'\n"
                "                       /|\\\n"
                "                        |\n"
                "               Strand I |\n"
                "                        |\n"
                "                        |\n"
                "                        |\n"
                "                    i1  5'\n\n");
    fprintf(fp, "     step      i1-i2\n");
  }

  k = (ds == 1) ? 1 : 4;
  for (i = 1; i < num_bp; i++) {
    init_dvector(olarea, 1, 10, 0.0);

    get_zoave(i, ds, orien, org, oave, zave);

    r1 = pair_num[1][i];     /* i1 */
    r2 = pair_num[1][i + 1]; /* i2 */

    if (with_stacking(ds, i, pair_num, orien, ring_atom, xyz)) {
      olarea[1] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 1);
      olarea[6] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 0);

      if (ds == 2) {
        r2 = pair_num[ds][i + 1]; /* j2 */
        olarea[2] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 1);
        olarea[7] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 0);
        r1 = pair_num[ds][i]; /* j1 */
        olarea[4] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 1);
        olarea[9] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 0);
        r2 = pair_num[1][i + 1]; /* i2 */
        olarea[3] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 1);
        olarea[8] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 0);
      }

      for (j = 1; j <= k; j++) {
        olarea[5] += olarea[j];
        olarea[10] += olarea[5 + j];
      }
    }

    if (ds == 1)
      fprintf(fp, "%4ld  %c/%c ", i, bp_seq[ds][i], bp_seq[ds][i + 1]);
    else
      fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
              bp_seq[ds][i + 1], bp_seq[ds][i]);
    for (j = 1; j <= k; j++)
      fprintf(fp, fmt, olarea[5 + j], olarea[j]);
    if (ds == 2)
      fprintf(fp, fmt, olarea[10], olarea[5]);
    fprintf(fp, "\n");
  }

  free_lmatrix(ring_atom, 1, num_residue, 1, 19);
  free_dvector(olarea, 1, 10);
}

/* Get xyz coordinates for base ring + 1st shell, offset by "oave" */
long ratom_xyz(long *ratom_list, long only_ring, double **xyz, double *oave,
               double **oxyz) {
  long i, k, n;

  n = ratom_list[10];
  for (i = 1; i <= n; i++) {
    k = (only_ring) ? ratom_list[i] : ratom_list[10 + i];
    ddxyz(oave, xyz[k], oxyz[i]);
  }

  return n;
}

/* Get average origin and z-axis for program "analyze" */
void get_zoave(long istep, long ds, double **orien, double **org, double *oave,
               double *zave) {
  double d;
  long ioffset3, ioffset9, j;

  ioffset3 = (istep - 1) * 3;
  ioffset9 = (istep - 1) * 9;

  d = dot(&orien[1][ioffset9 + 6], &orien[1][ioffset9 + 15]); /* i1 vs i2 */
  avexyz(org[1] + ioffset3, org[1] + ioffset3 + 3, oave);
  (d > 0.0) ? sumxyz(orien[1] + ioffset9 + 15, orien[1] + ioffset9 + 6, zave)
            : ddxyz(orien[1] + ioffset9 + 15, orien[1] + ioffset9 + 6, zave);

  if (ds == 2) {
    for (j = 1; j <= 3; j++) {
      oave[j] += 0.5 * (org[ds][ioffset3 + j] + org[ds][ioffset3 + 3 + j]);
      oave[j] *= 0.5;
    }
    d = dot(&orien[1][ioffset9 + 6], &orien[ds][ioffset9 + 6]); /* i1 vs j1 */
    for (j = 1; j <= 3; j++)
      zave[j] += (d > 0.0) ? orien[ds][ioffset9 + 6 + j]
                           : -orien[ds][ioffset9 + 6 + j];
    d = dot(&orien[1][ioffset9 + 6], &orien[ds][ioffset9 + 15]); /* i1 vs j2 */
    for (j = 1; j <= 3; j++)
      zave[j] += (d > 0.0) ? orien[ds][ioffset9 + 15 + j]
                           : -orien[ds][ioffset9 + 15 + j];
  }
  vec_norm(zave);
}

/* Get average origin and z-axis for a base-pair */
void get_bp_zoave(long ia, long ib, double **orien, double **org, double *oave,
                  double *zave) {
  double d;

  d = dot(&orien[ia][6], &orien[ib][6]);
  avexyz(org[ia], org[ib], oave);
  (d > 0.0) ? sumxyz(orien[ib] + 6, orien[ia] + 6, zave)
            : ddxyz(orien[ib] + 6, orien[ia] + 6, zave);
  vec_norm(zave);
}

/* Get ring atom index and linkage information in base residues: exclude Hs */
void ring_oidx(long num, long num_residue, long *RY, long **seidx,
               char **AtomName, double **xyz, long *idx, long **ring_atom) {
  long i, num_ring;
  long **connect;

  connect = lmatrix(1, num, 1, 7); /* overall connection */
  get_bonds(num, AtomName, xyz, num_residue, RY, seidx, connect);
  all_bring_atoms(num_residue, RY, seidx, AtomName, &num_ring, ring_atom);

  for (i = 1; i <= num_residue; i++) /* for each residue */
    if (ring_atom[i][10] > 0)
      get_cntatom(ring_atom[i], connect, idx);

  free_lmatrix(connect, 1, num, 1, 7);
}

/* Get ring atom connections */
void get_cntatom(long *ringlist, long **connect, long *idx) {
  long i, ic = 0, id, ix, j, ra_num;

  ra_num = ringlist[10]; /* # of ring atoms */
  for (i = 1; i <= ra_num; i++) {
    id = ringlist[i];
    ix = 0;
    for (j = 1; j <= connect[id][7]; j++) {
      ic = connect[id][j];
      if (idx[ic] == 3)
        continue;                                  /* excluding H atom */
      if (!lval_in_set(ic, 1, ra_num, ringlist)) { /* not a ring atom */
        ix = 1;
        break;
      }
    }
    ringlist[10 + i] = (ix) ? ic : id;
  }
}

/* ============ calculating polygon overlap area ============ */
static void pia_fit(double minx, double miny, const double mid, double sclx,
                    double scly, point *x, long cx, vertex *ix, long fudge) {
  long c, t;
  point t1, t2;

  c = cx;
  while (c--) {
    t = (x[c].x - minx) * sclx - mid;
    ix[c].ip.x = (t & ~7) | fudge | (c & 1);
    t = (x[c].y - miny) * scly - mid;
    ix[c].ip.y = (t & ~7) | fudge;
  }

  ix[0].ip.y += cx & 1;
  ix[cx] = ix[0];

  c = cx;
  while (c--) {
    t1.x = ix[c].ip.x;
    t1.y = ix[c + 1].ip.x;
    t2.x = ix[c + 1].ip.x;
    t2.y = ix[c].ip.x;
    ix[c].rx = (ix[c].ip.x < ix[c + 1].ip.x) ? t1 : t2;

    t1.x = ix[c].ip.y;
    t1.y = ix[c + 1].ip.y;
    t2.x = ix[c + 1].ip.y;
    t2.y = ix[c].ip.y;
    ix[c].ry = (ix[c].ip.y < ix[c + 1].ip.y) ? t1 : t2;
    ix[c].inside = 0;
  }
}

static double pia_area(point a, point p, point q) {
  return p.x * q.y - p.y * q.x + a.x * (p.y - q.y) + a.y * (q.x - p.x);
}

static void pia_cntrib(double *s, point f, point t, long w) {
  *s += w * (t.x - f.x) * (t.y + f.y) * 0.5;
}

static long pia_ovl(point p, point q) { return p.x < q.y && q.x < p.y; }

static void pia_cross(double *out_s, vertex *a, vertex *b, vertex *c, vertex *d,
                      double a1, double a2, double a3, double a4) {
  double r1, r2;
  point dp;

  r1 = a1 / (a1 + a2);
  r2 = a3 / (a3 + a4);

  dp.x = a->ip.x + r1 * (b->ip.x - a->ip.x);
  dp.y = a->ip.y + r1 * (b->ip.y - a->ip.y);
  pia_cntrib(out_s, dp, b->ip, 1);

  dp.x = c->ip.x + r2 * (d->ip.x - c->ip.x);
  dp.y = c->ip.y + r2 * (d->ip.y - c->ip.y);
  pia_cntrib(out_s, d->ip, dp, 1);

  ++a->inside;
  --c->inside;
}

static void pia_inness(double *out_s, vertex *P, long cP, vertex *Q, long cQ) {
  long j, sgn, s = 0, c = cQ;
  point p = P[0].ip;

  while (c--)
    if (Q[c].rx.x < p.x && p.x < Q[c].rx.y) {
      sgn = 0 < pia_area(p, Q[c].ip, Q[c + 1].ip);
      s += sgn != (Q[c].ip.x < Q[c + 1].ip.x) ? 0 : (sgn ? -1 : 1);
    }
  for (j = 0; j < cP; ++j) {
    if (s)
      pia_cntrib(out_s, P[j].ip, P[j + 1].ip, s);
    s += P[j].inside;
  }
}

/* Area of intersection between two polygons */
static double pia_inter(point *a, long na, point *b, long nb) {
  double minx, miny, maxx, maxy, ascale, sclx, scly;
  double out_s = 0.0, a1, a2, a3, a4;
  const double gamut = 5.0e8, mid = 0.5 * gamut;
  long j, k, o;
  vertex ipa[MNPOLY + 1],
      ipb[MNPOLY + 1]; /* ISO C89 forbids variable-size array */

  if (na < 3 || nb < 3)
    return 0;

  minx = miny = XBIG;
  maxx = maxy = -XBIG;
  for (j = 0; j < na; j++) {
    if (minx > a[j].x)
      minx = a[j].x;
    if (miny > a[j].y)
      miny = a[j].y;
    if (maxx < a[j].x)
      maxx = a[j].x;
    if (maxy < a[j].y)
      maxy = a[j].y;
  }
  for (j = 0; j < nb; j++) {
    if (minx > b[j].x)
      minx = b[j].x;
    if (miny > b[j].y)
      miny = b[j].y;
    if (maxx < b[j].x)
      maxx = b[j].x;
    if (maxy < b[j].y)
      maxy = b[j].y;
  }

  sclx = gamut / (maxx - minx);
  scly = gamut / (maxy - miny);
  ascale = sclx * scly;

  pia_fit(minx, miny, mid, sclx, scly, a, na, ipa, 0);
  pia_fit(minx, miny, mid, sclx, scly, b, nb, ipb, 2);

  for (j = 0; j < na; ++j)
    for (k = 0; k < nb; ++k)
      if (pia_ovl(ipa[j].rx, ipb[k].rx) && pia_ovl(ipa[j].ry, ipb[k].ry)) {
        a1 = -pia_area(ipa[j].ip, ipb[k].ip, ipb[k + 1].ip);
        a2 = pia_area(ipa[j + 1].ip, ipb[k].ip, ipb[k + 1].ip);
        o = a1 < 0;
        if (o == (a2 < 0)) {
          a3 = pia_area(ipb[k].ip, ipa[j].ip, ipa[j + 1].ip);
          a4 = -pia_area(ipb[k + 1].ip, ipa[j].ip, ipa[j + 1].ip);
          if ((a3 < 0) == (a4 < 0)) {
            if (o)
              pia_cross(&out_s, &ipa[j], &ipa[j + 1], &ipb[k], &ipb[k + 1], a1,
                        a2, a3, a4);
            else
              pia_cross(&out_s, &ipb[k], &ipb[k + 1], &ipa[j], &ipa[j + 1], a3,
                        a4, a1, a2);
          }
        }
      }
  pia_inness(&out_s, ipa, na, ipb, nb);
  pia_inness(&out_s, ipb, nb, ipa, na);

  return fabs(out_s) / ascale;
}

/* ============ calculating polygon overlap area ============ */

/* Get overlap area between base residues r1 and r2 */
double get_oarea(long r1, long r2, long **ring_atom, double *oave, double *zave,
                 double **xyz, long only_ring) {
  long i, j, n1, n2;
  double **oxyz1, **oxyz2, **oxyz1Z, **oxyz2Z, **rotmat;
  point a[MNPOLY], b[MNPOLY];

  oxyz1 = dmatrix(1, 9, 1, 3);
  oxyz2 = dmatrix(1, 9, 1, 3);
  oxyz1Z = dmatrix(1, 9, 1, 3);
  oxyz2Z = dmatrix(1, 9, 1, 3);
  rotmat = dmatrix(1, 3, 1, 3); /* no use here */

  n1 = ratom_xyz(ring_atom[r1], only_ring, xyz, oave, oxyz1);
  n2 = ratom_xyz(ring_atom[r2], only_ring, xyz, oave, oxyz2);

  align2zaxis(n1, zave, rotmat, oxyz1, oxyz1Z);
  align2zaxis(n2, zave, rotmat, oxyz2, oxyz2Z);

  /* change xy1 & xy2 to an array of point */
  for (i = 1; i <= n1; i++) {
    j = i - 1;
    a[j].x = oxyz1Z[i][1];
    a[j].y = oxyz1Z[i][2];
  }
  for (i = 1; i <= n2; i++) {
    j = i - 1;
    b[j].x = oxyz2Z[i][1];
    b[j].y = oxyz2Z[i][2];
  }

  free_dmatrix(oxyz1, 1, 9, 1, 3);
  free_dmatrix(oxyz2, 1, 9, 1, 3);
  free_dmatrix(oxyz1Z, 1, 9, 1, 3);
  free_dmatrix(oxyz2Z, 1, 9, 1, 3);
  free_dmatrix(rotmat, 1, 3, 1, 3);

  return pia_inter(a, n1, b, n2);
}

/* Calculate the six local CEHS base-pair parameters.
   Propeller is applied first followed by buckle-opening */
void cehs_bppar(double **rot1, double *org1, double **rot2, double *org2,
                double *pars, double **mst_orien, double *mst_org) {
  double buckleopening, phi;
  double hinge[4], t1[4], t2[4], xm[4], ym[4], zm[4];
  double **paraII, **paraI, **temp;
  long i, j;

  for (i = 1; i <= 3; i++) {
    t1[i] = rot1[i][2]; /* y1 */
    t2[i] = rot2[i][2]; /* y2 */
  }

  cross(t1, t2, hinge);
  buckleopening = magang(t1, t2);

  temp = dmatrix(1, 3, 1, 3);
  paraII = dmatrix(1, 3, 1, 3);
  paraI = dmatrix(1, 3, 1, 3);

  arb_rotation(hinge, -0.5 * buckleopening, temp);
  multi_matrix(temp, 3, 3, rot2, 3, 3, paraI);
  arb_rotation(hinge, 0.5 * buckleopening, temp);
  multi_matrix(temp, 3, 3, rot1, 3, 3, paraII);

  for (i = 1; i <= 3; i++) {
    ym[i] = paraI[i][2];  /* also paraII[i][2] */
    t1[i] = paraII[i][1]; /* x1 */
    t2[i] = paraI[i][1];  /* x2 */
  }

  /* twist is the angle between the two y- or x-axes */
  pars[5] = vec_ang(t1, t2, ym);

  sumxyz(t1, t2, xm);
  vec_norm(xm);

  cross(xm, ym, zm);

  avexyz(org1, org2, mst_org);
  ddxyz(org1, org2, t1);
  x_y_z_2_mtx(xm, ym, zm, mst_orien);

  /* get the xyz displacement parameters */
  for (i = 1; i <= 3; i++) {
    pars[i] = 0.0;
    for (j = 1; j <= 3; j++)
      pars[i] += t1[j] * mst_orien[j][i];
  }

  /* phi angle is defined by hinge and xm */
  phi = deg2rad(vec_ang(hinge, xm, ym));

  /* get buckle and opening angles */
  pars[4] = buckleopening * cos(phi);
  pars[6] = buckleopening * sin(phi);

  free_dmatrix(temp, 1, 3, 1, 3);
  free_dmatrix(paraII, 1, 3, 1, 3);
  free_dmatrix(paraI, 1, 3, 1, 3);
}

static void check_cehs_pair_frames(double **r1, double **r2) {
  double sum = 0;
  long i;

  /* sum = dot(z1, z2) */
  for (i = 1; i <= 3; i++)
    sum += r1[i][3] * r2[i][3];

  if (sum > 0)
    return;

  /* opposite direction, as in a WC pair */
  for (i = 1; i <= 3; i++) {
    r1[i][2] = -r1[i][2]; /* reverse y */
    r1[i][3] = -r1[i][3]; /* reverse z */
  }
}

/* In CEHS: for base-pair parameters, propeller is applied first */
void out_cehs(long num_bp, char **bp_seq, long *bphlx, double **orien,
              double **org, FILE *fp) {
  double mfoi[4], o1[4], o2[4];
  double *bp_org, *bp_orien, **bp_par, **mfi, **r1, **r2, **step_par;
  long i, nbpm1;
  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);

  bp_org = dvector(1, num_bp * 3);
  bp_orien = dvector(1, num_bp * 9);

  bp_par = dmatrix(1, num_bp, 1, 6);
  step_par = dmatrix(1, nbpm1, 1, 6);

  for (i = 1; i <= num_bp; i++) {
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    cehs_bppar(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

    cpxyz(mfoi, bp_org + (i - 1) * 3);
    mst2orien(bp_orien, (i - 1) * 9, mfi);
  }

  for (i = 1; i <= nbpm1; i++) {
    refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
    bpstep_par(r1, o1, r2, o2, step_par[i], mfi, mfoi);
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "CEHS base-pair parameters\n");
  print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

  print_sep(fp, '*', 76);
  fprintf(fp, "CEHS base-pair step parameters\n");
  prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dvector(bp_org, 1, num_bp * 3);
  free_dvector(bp_orien, 1, num_bp * 9);
  free_dmatrix(bp_par, 1, num_bp, 1, 6);
  free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters based on Gorin's scheme */
void compdna(double **rot1, double *org1, double **rot2, double *org2,
             double *pars, double **mst_orien, double *mst_org) {
  double dorg[4], dz[4], y1_osx[4], y2[4], z1[4], z2[4];
  double xm[4], ym[4], zm[4];
  long i;

  /* use only z- and y-axes for constructing middle frame */
  for (i = 1; i <= 3; i++) {
    z1[i] = rot1[i][3];
    y1_osx[i] = rot1[i][2];
    z2[i] = rot2[i][3];
    y2[i] = rot2[i][2];
  }
  avexyz(org1, org2, mst_org);
  ddxyz(org1, org2, dorg);
  sumxyz(z1, z2, zm);
  ddxyz(z1, z2, dz);
  vec_norm(zm);

  vec_orth(y1_osx, zm); /* orthogonal y-component */
  vec_orth(y2, zm);
  sumxyz(y1_osx, y2, ym);
  vec_norm(ym);

  cross(ym, zm, xm);

  pars[1] = dot(dorg, xm);
  pars[2] = dot(dorg, ym);
  pars[3] = dot(dorg, zm);
  pars[4] = -2 * rad2deg(asin(dot(dz, ym) / 2));
  pars[5] = 2 * rad2deg(asin(dot(dz, xm) / 2));
  pars[6] = vec_ang(y1_osx, y2, zm);

  x_y_z_2_mtx(xm, ym, zm, mst_orien);
}

void out_compdna(long num_bp, char **bp_seq, long *bphlx, double **orien,
                 double **org, FILE *fp) {
  double mfoi[4], o1[4], o2[4];
  double *bp_org, *bp_orien, **bp_par, **mfi, **r1, **r2, **step_par;
  long i, nbpm1;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);

  bp_org = dvector(1, num_bp * 3);
  bp_orien = dvector(1, num_bp * 9);

  bp_par = dmatrix(1, num_bp, 1, 6);
  step_par = dmatrix(1, nbpm1, 1, 6);

  for (i = 1; i <= num_bp; i++) {
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    compdna(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

    cpxyz(mfoi, bp_org + (i - 1) * 3);
    mst2orien(bp_orien, (i - 1) * 9, mfi);
  }

  for (i = 1; i <= nbpm1; i++) {
    refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
    compdna(r1, o1, r2, o2, step_par[i], mfi, mfoi);
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "CompDNA base-pair parameters\n");
  print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

  print_sep(fp, '*', 76);
  fprintf(fp, "CompDNA base-pair step parameters\n");
  prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dvector(bp_org, 1, num_bp * 3);
  free_dvector(bp_orien, 1, num_bp * 9);
  free_dmatrix(bp_par, 1, num_bp, 1, 6);
  free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Calculate Curves parameters */
void my_curves(double **rot1, double *org1, double **rot2, double *org2,
               double *pars) {
  double dl, du, twist_l, twist_u;
  double j1_osx[4], j2[4], k1[4], k2[4], l1[4], l2[4];
  double d[4], f[4], n[4], pl[4], pu[4], q[4];
  double temp[4], rtmp[4];
  long i;

  /* decompose rot1 and rot2 into JKL */
  for (i = 1; i <= 3; i++) {
    j1_osx[i] = rot1[i][1];
    k1[i] = rot1[i][2];
    l1[i] = rot1[i][3];
    j2[i] = rot2[i][1];
    k2[i] = rot2[i][2];
    l2[i] = rot2[i][3];
  }
  avexyz(org1, org2, q); /* mean origin */
  sumxyz(l1, l2, n);     /* n is z-axis */
  sumxyz(j1_osx, j2, d); /* d is x-axis */
  vec_norm(n);
  vec_orth(d, n);
  cross(n, d, f); /* f is y-axis */

  /* get the intersection of l1 & l2 with the above mean-plane */
  ddxyz(org1, q, temp);
  dl = dot(n, temp) / dot(n, l1); /* org1 to mean-plane along l1 */
  for (i = 1; i <= 3; i++)
    pl[i] = org1[i] + dl * l1[i]; /* l1 intersection with mean-plane */

  ddxyz(q, org2, temp);
  du = dot(n, temp) / dot(n, l2); /* org2 to mean-plane along l2 */
  for (i = 1; i <= 3; i++)
    pu[i] = org2[i] - du * l2[i]; /* l2 intersection with mean-plane */

  pars[3] = dl + du;      /* this is why RISE is bigger */
  ddxyz(pl, pu, temp);    /* vector pl ---> pu */
  pars[1] = dot(temp, d); /* Shift */
  pars[2] = dot(temp, f); /* Slide */

  cross(l2, d, temp);
  pars[4] = 2 * vec_ang(f, temp, d); /* Tilt */

  cross(d, temp, rtmp);
  pars[5] = 2 * vec_ang(rtmp, l2, temp); /* Roll */

  get_vector(f, d, -pars[4] / 2, temp);
  twist_l = vec_ang(k1, temp, l1);

  get_vector(f, d, +pars[4] / 2, temp);
  twist_u = vec_ang(temp, k2, l2);

  pars[6] = twist_l + twist_u; /* Twist */
}

/* Get mean base-pair frame using Curves method */
void curves_mbt(long ibp, double **orien, double **org, double **cvr,
                double *cvo) {
  double o1[4], o2[4], xm[4], ym[4], zm[4];
  double **r1, **r2;
  long i;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);

  refs_right_left(ibp, orien, org, r1, o1, r2, o2);
  for (i = 1; i <= 3; i++) {
    zm[i] = r1[i][3] + r2[i][3];
    xm[i] = r1[i][1] + r2[i][1];
  }
  vec_norm(zm);
  vec_norm(xm);
  cross(zm, xm, ym); /* xm & zm not orthogonal */
  x_y_z_2_mtx(xm, ym, zm, cvr);
  avexyz(o1, o2, cvo);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
}

void out_curves(long num_bp, char **bp_seq, long *bphlx, double **orien,
                double **org, FILE *fp) {
  double o1[4], o2[4];
  double **bp_par, **mfi, **r1, **r2, **step_par;
  long i, nbpm1;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);

  bp_par = dmatrix(1, num_bp, 1, 6);
  step_par = dmatrix(1, nbpm1, 1, 6);

  for (i = 1; i <= num_bp; i++) {
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    my_curves(r1, o1, r2, o2, bp_par[i]);
  }

  for (i = 1; i <= nbpm1; i++) {
    curves_mbt(i, orien, org, r1, o1);
    curves_mbt(i + 1, orien, org, r2, o2);
    my_curves(r1, o1, r2, o2, step_par[i]);
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "Curves base-pair parameters (by Xiang-Jun Lu)\n");
  print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

  print_sep(fp, '*', 76);
  fprintf(fp, "Curves base-pair step parameters\n");
  prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dmatrix(bp_par, 1, num_bp, 1, 6);
  free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters based on Dickerson's scheme */
void freehelix(double **rot1, double *org1, double **rot2, double *org2,
               double *pars, double **mst_orien, double *mst_org) {
  double dorg[4], y1_osx[4], y2[4], z1[4], z2[4];
  double xm[4], ym[4], zm[4];
  long i;

  /* use only y- and z-axes for constructing middle frame */
  for (i = 1; i <= 3; i++) {
    y1_osx[i] = rot1[i][2];
    z1[i] = rot1[i][3];
    y2[i] = rot2[i][2];
    z2[i] = rot2[i][3];
  }
  avexyz(org1, org2, mst_org);
  ddxyz(org1, org2, dorg);
  sumxyz(y1_osx, y2, ym);
  sumxyz(z1, z2, zm);
  vec_norm(ym);

  vec_norm(zm); /* mst-z */
  cross(ym, zm, xm);
  vec_norm(xm);      /* mst-x */
  cross(zm, xm, ym); /* mst-y */
  vec_norm(ym);

  pars[1] = dot(dorg, xm);
  pars[2] = dot(dorg, ym);
  pars[3] = dot(dorg, zm);
  pars[4] = vec_ang(z1, z2, xm);
  pars[5] = vec_ang(z1, z2, ym);
  pars[6] = vec_ang(y1_osx, y2, zm);

  x_y_z_2_mtx(xm, ym, zm, mst_orien);
}

void out_freehelix(long num_bp, char **bp_seq, long *bphlx, double **orien,
                   double **org, FILE *fp) {
  double mfoi[4], o1[4], o2[4];
  double *bp_org, *bp_orien, **bp_par, **mfi, **r1, **r2, **step_par;
  long i, nbpm1;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);

  bp_org = dvector(1, num_bp * 3);
  bp_orien = dvector(1, num_bp * 9);

  bp_par = dmatrix(1, num_bp, 1, 6);
  step_par = dmatrix(1, nbpm1, 1, 6);

  for (i = 1; i <= num_bp; i++) {
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    freehelix(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

    cpxyz(mfoi, bp_org + (i - 1) * 3);
    mst2orien(bp_orien, (i - 1) * 9, mfi);
  }

  for (i = 1; i <= nbpm1; i++) {
    refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
    freehelix(r1, o1, r2, o2, step_par[i], mfi, mfoi);
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "FreeHelix base-pair parameters (by Xiang-Jun Lu)\n");
  print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

  print_sep(fp, '*', 76);
  fprintf(fp, "FreeHelix base-pair step parameters\n");
  prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dvector(bp_org, 1, num_bp * 3);
  free_dvector(bp_orien, 1, num_bp * 9);
  free_dmatrix(bp_par, 1, num_bp, 1, 6);
  free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Get the single helical rotation angle and its corresponding axis */
void sgl_helix(double **rot1, double **rot2, double *rot_ang, double *rot_hlx) {
  double dsum = 0.0, tcos;
  double dx[4] = {EMPTY_NUMBER, 1.0, 0.0, 0.0};
  double dy[4] = {EMPTY_NUMBER, 0.0, 1.0, 0.0};
  double **R, **temp;
  long i, ichg = 0, j;

  R = dmatrix(1, 3, 1, 3);
  temp = dmatrix(1, 3, 1, 3);

  transpose_matrix(rot1, 3, 3, temp);
  multi_matrix(temp, 3, 3, rot2, 3, 3, R); /* rotation matrix w.r.t. 1 */

  for (i = 1; i <= 3; i++)
    dsum += R[i][i];         /* trace of R */
  tcos = 0.5 * (dsum - 1.0); /* positive rotation angle */
  *rot_ang = (tcos >= 1.0) ? 0.0 : rad2deg(acos(tcos));

  /* helical rotation axis, also = cross(x1-x2, y1-y2) = cross(x2-x1, y2-y1) */
  for (i = 1; i <= 3; i++) {
    dx[i] = R[i][1] - dx[i];
    dy[i] = R[i][2] - dy[i];
  }
  cross(dx, dy, rot_hlx);
  vec_norm(rot_hlx);

  /* check back */
  arb_rotation(rot_hlx, *rot_ang, temp);
  for (i = 1; i <= 3 && !ichg; i++) {
    for (j = 1; j <= 3 && !ichg; j++)
      if (fabs(R[i][j] - temp[i][j]) > XEPS) {
        *rot_ang = -*rot_ang; /* negate rotation angle */
        ichg = 1;
      }
  }

#if 0 /* September 22, 2006: for verification */
{
        if (ichg) {
            for (i = 1; i <= 3; i++) {
                for (j = 1; j <= 3; j++)
                    fprintf(stderr, "\t%8.5f(%8.5f)", R[i][j], temp[i][j]);
                fprintf(stderr, "\n");
            }
            print_sep(stderr, '-', 78);
            arb_rotation(rot_hlx, *rot_ang, temp);
            for (i = 1; i <= 3; i++) {
                for (j = 1; j <= 3; j++)
                    fprintf(stderr, "\t%8.5f(%8.5f)", R[i][j], temp[i][j]);
                fprintf(stderr, "\n");
            }
        }
    }
#endif

  free_dmatrix(R, 1, 3, 1, 3);
  free_dmatrix(temp, 1, 3, 1, 3);
}

/* Calculate parameters based on Tung's scheme */
void ngeom(double **rot1, double *org1, double **rot2, double *org2,
           double *pars, double **mst_orien, double *mst_org) {
  double ang;
  double dorg[4], dorg1[4], haxis[4], tpars[7];
  double **rmtx;
  long i;

  rmtx = dmatrix(1, 3, 1, 3);

  /* translational parameters are the same as in RNA */
  sgl_helix(rot1, rot2, &ang, haxis);

  ddxyz(org1, org2, dorg);
  multi_vec_matrix(dorg, 3, rot1, 3, 3, dorg1);
  arb_rotation(haxis, ang / 2, rmtx);
  multi_vec_matrix(dorg1, 3, rmtx, 3, 3, pars);

  /* rotational parameters are the same as in CEHS */
  bpstep_par(rot1, org1, rot2, org2, tpars, mst_orien, mst_org);

  for (i = 4; i <= 6; i++)
    pars[i] = tpars[i];

  free_dmatrix(rmtx, 1, 3, 1, 3);
}

void out_ngeom(long num_bp, char **bp_seq, long *bphlx, double **orien,
               double **org, FILE *fp) {
  double mfoi[4], o1[4], o2[4];
  double *bp_org, *bp_orien, **bp_par, **mfi, **r1, **r2, **step_par;
  long i, nbpm1;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);

  bp_org = dvector(1, num_bp * 3);
  bp_orien = dvector(1, num_bp * 9);

  bp_par = dmatrix(1, num_bp, 1, 6);
  step_par = dmatrix(1, nbpm1, 1, 6);
  for (i = 1; i <= num_bp; i++) {
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    ngeom(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

    cpxyz(mfoi, bp_org + (i - 1) * 3);
    mst2orien(bp_orien, (i - 1) * 9, mfi);
  }

  for (i = 1; i <= nbpm1; i++) {
    refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
    ngeom(r1, o1, r2, o2, step_par[i], mfi, mfoi);
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "NGEOM base-pair parameters\n");
  print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

  print_sep(fp, '*', 76);
  fprintf(fp, "NGEOM base-pair step parameters\n");
  prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dvector(bp_org, 1, num_bp * 3);
  free_dvector(bp_orien, 1, num_bp * 9);
  free_dmatrix(bp_par, 1, num_bp, 1, 6);
  free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters based on Bansal's scheme */
void nuparm(double **rot1, double *org1, double **rot2, double *org2,
            double *pars, double **mst_orien, double *mst_org, double *hpars,
            long get_hpar) {
  double a, cx, cy, cz, dx, dy, sx, sy, sz;
  double B[4], dorg[4], x1[4], x2[4], y1_osx[4], y2[4], zh[4];
  double xm[4], ym[4], zm[4];
  double **A, **invA;
  long i;

  for (i = 1; i <= 3; i++) {
    x1[i] = rot1[i][1];
    y1_osx[i] = rot1[i][2];
    x2[i] = rot2[i][1];
    y2[i] = rot2[i][2];
  }
  sumxyz(x1, x2, xm);
  sumxyz(y1_osx, y2, ym);
  avexyz(org1, org2, mst_org);
  ddxyz(org1, org2, dorg);

  /* get the middle frame (xm & ym are not orthogonal) */
  vec_norm(xm);
  vec_norm(ym);
  cross(xm, ym, zm);
  vec_norm(zm);

  x_y_z_2_mtx(xm, ym, zm, mst_orien);

  pars[1] = dot(dorg, xm);
  pars[2] = dot(dorg, ym);
  pars[3] = dot(dorg, zm);
  pars[4] = -2 * rad2deg(asin(dot(y1_osx, zm)));
  pars[5] = 2 * rad2deg(asin(dot(x1, zm)));
  pars[6] = vec_ang(y1_osx, y2, zm);

  if (get_hpar) { /* helical parameters */
    ddxyz(x2, x1, xm);
    ddxyz(y2, y1_osx, ym);
    cross(xm, ym, zh);
    vec_norm(zh);

    hpars[6] = vec_ang(y1_osx, y2, zh);
    hpars[5] = -rad2deg(asin(dot(zh, x1)));
    hpars[4] = rad2deg(asin(dot(zh, y1_osx)));

    a = deg2rad(hpars[4]);
    cx = cos(a);
    sx = sin(a);
    a = deg2rad(hpars[5]);
    cy = cos(a);
    sy = -sin(a);
    a = deg2rad(pars[6]); /* not helical twist */
    cz = cos(a);
    sz = sin(a);

    A = dmatrix(1, 3, 1, 3);
    invA = dmatrix(1, 3, 1, 3);

    A[1][1] = 2 * cx * sz;
    A[1][2] = 0.0;
    A[1][3] = 2 * sx;
    A[2][1] = 0.0;
    A[2][2] = -2 * cy * sz;
    A[2][3] = 2 * sy;
    A[3][1] = -2 * cy * sx * sz;
    A[3][2] = 2 * cx * sy * sz;
    A[3][3] = cx * cy * (1 + cz);
    dinverse(A, 3, invA);
    transpose_matrix(invA, 3, 3, A);

    dx = sqrt(2 * (1 + cz + sx * sx * (1 - cz)));
    dy = sqrt(2 * (1 + cz + sy * sy * (1 - cz)));
    B[1] = pars[2] * dy;
    B[2] = pars[1] * dx;
    B[3] = 0.5 * pars[3] * dx * dx;

    multi_vec_matrix(B, 3, A, 3, 3, dorg);
    cpxyz(dorg, hpars);

    free_dmatrix(A, 1, 3, 1, 3);
    free_dmatrix(invA, 1, 3, 1, 3);
  }
}

void out_nuparm(long num_bp, char **bp_seq, long *bphlx, double **orien,
                double **org, FILE *fp) {
  double hpi[7], mfoi[4], o1[4], o2[4];
  double *bp_org, *bp_orien, **bp_par, **mfi, **heli_par, **r1, **r2,
      **step_par;
  long i, nbpm1;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);

  bp_org = dvector(1, num_bp * 3);
  bp_orien = dvector(1, num_bp * 9);

  bp_par = dmatrix(1, num_bp, 1, 6);
  step_par = dmatrix(1, nbpm1, 1, 6);
  heli_par = dmatrix(1, nbpm1, 1, 6);

  for (i = 1; i <= num_bp; i++) {
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    nuparm(r1, o1, r2, o2, bp_par[i], mfi, mfoi, hpi, 0);

    cpxyz(mfoi, bp_org + (i - 1) * 3);
    mst2orien(bp_orien, (i - 1) * 9, mfi);
  }

  for (i = 1; i <= nbpm1; i++) {
    refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
    nuparm(r1, o1, r2, o2, step_par[i], mfi, mfoi, heli_par[i], 1);
  }

  print_sep(fp, '*', 76);
  fprintf(fp, "NUPARM base-pair parameters\n");
  print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

  print_sep(fp, '*', 76);
  fprintf(fp, "NUPARM base-pair step parameters\n");
  prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

  print_sep(fp, '*', 76);
  fprintf(fp, "NUPARM base-pair helical parameters\n");
  prt_step_par(bp_seq, num_bp, bphlx, 1, heli_par, fp);

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dvector(bp_org, 1, num_bp * 3);
  free_dvector(bp_orien, 1, num_bp * 9);
  free_dmatrix(bp_par, 1, num_bp, 1, 6);
  free_dmatrix(step_par, 1, nbpm1, 1, 6);
  free_dmatrix(heli_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters based on Babcock's scheme */
void rna(double **rot1, double *org1, double *pvt1, double **rot2, double *org2,
         double *pvt2, double *pars, double **mst_orien, double *mst_org) {
  double ang;
  double dorg[4], dorg1[4], haxis[4], p1[4], p2[4], pt[4];
  double **rmtx;
  long i;

  rmtx = dmatrix(1, 3, 1, 3);

  /* total rotation angle and helical axis */
  sgl_helix(rot1, rot2, &ang, haxis);

  /* rotational parameters */
  for (i = 1; i <= 3; i++)
    pars[i + 3] = ang * haxis[i];
  ddxyz(org1, org2, dorg);

  /* translational parameters */
  multi_vec_matrix(dorg, 3, rot1, 3, 3, dorg1);
  arb_rotation(haxis, ang / 2, rmtx);
  multi_matrix(rot1, 3, 3, rmtx, 3, 3, mst_orien);
  multi_vec_Tmatrix(pvt1, 3, rot1, 3, 3, p1);
  multi_vec_Tmatrix(pvt2, 3, rot2, 3, 3, p2);
  for (i = 1; i <= 3; i++)
    mst_org[i] = 0.5 * (org1[i] + org2[i] + p1[i] + p2[i]);
  ddxyz(pvt1, dorg1, pt);
  multi_vec_matrix(pt, 3, rmtx, 3, 3, p1);
  multi_vec_Tmatrix(pvt2, 3, rmtx, 3, 3, p2);
  for (i = 1; i <= 3; i++)
    pars[i] = p1[i] + p2[i] + pvt1[i] - pvt2[i];

  free_dmatrix(rmtx, 1, 3, 1, 3);
}

/* Make correction to dx and dy using Babcock's pivot point */
void pvt_dxdy(double **rot1, double *org1, double *pvt1, double *pars,
              double **mst_orien, double *mst_org) {
  double TipInc1, half_rise;
  double axis_h[4], hinge1[4], org1_h[4], t1[4], t2[4];
  double **rot1_h, **temp;
  long i, j;

  temp = dmatrix(1, 3, 1, 3);
  rot1_h = dmatrix(1, 3, 1, 3);

  half_rise = 0.5 * pars[3];
  for (i = 1; i <= 3; i++) {
    t1[i] = rot1[i][3];          /* z1 */
    axis_h[i] = mst_orien[i][3]; /* helical axis */
    org1_h[i] = mst_org[i] - half_rise * mst_orien[i][3];
  }
  TipInc1 = magang(axis_h, t1);
  cross(axis_h, t1, hinge1);
  arb_rotation(hinge1, -TipInc1, temp);
  multi_matrix(temp, 3, 3, rot1, 3, 3, rot1_h);

  for (i = 1; i <= 3; i++) {
    t1[i] = dot(pvt1, rot1[i]);
    negate_xyz(temp[i]);
    temp[i][i] += 1.0;
  }
  for (i = 1; i <= 3; i++)
    t2[i] = dot(t1, temp[i]) + org1[i] - org1_h[i];

  for (i = 1; i <= 2; i++) {
    pars[i] = 0.0;
    for (j = 1; j <= 3; j++)
      pars[i] += t2[j] * rot1_h[j][i];
  }

  free_dmatrix(rot1_h, 1, 3, 1, 3);
  free_dmatrix(temp, 1, 3, 1, 3);
}

void out_rna(long ds, long num_bp, char **bp_seq, long *bphlx, double **orien,
             double **org, FILE *fp) {
  double pvt0[4] = {EMPTY_NUMBER, 0.0, 0.0, 0.0};
  double pvt1[4] = {EMPTY_NUMBER, 0.0, 1.808, 0.0};
  double pvt2[4] = {EMPTY_NUMBER, 0.0, -1.808, 0.0};
  double hpi[7], mfoi[4], o1[4], o2[4];
  double **mfi, **r1, **r2, **heli_par, **step_par;
  long i, j, k, nbpm1;

  nbpm1 = num_bp - 1;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mfi = dmatrix(1, 3, 1, 3);

  step_par = dmatrix(1, nbpm1, 1, 6);
  heli_par = dmatrix(1, nbpm1, 1, 6);

  if (ds == 2) { /* double helix */
    double *bp_org, *bp_orien, **bp_par;

    bp_org = dvector(1, num_bp * 3);
    bp_orien = dvector(1, num_bp * 9);
    bp_par = dmatrix(1, num_bp, 1, 6);

    for (i = 1; i <= num_bp; i++) {
      refs_right_left(i, orien, org, r1, o1, r2, o2);
      rna(r1, o1, pvt2, r2, o2, pvt1, bp_par[i], mfi, mfoi);

      cpxyz(mfoi, bp_org + (i - 1) * 3);
      mst2orien(bp_orien, (i - 1) * 9, mfi);
    }

    for (i = 1; i <= nbpm1; i++) {
      refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
      rna(r1, o1, pvt0, r2, o2, pvt0, step_par[i], mfi, mfoi);
      helical_par(r1, o1, r2, o2, heli_par[i], mfi, mfoi);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "RNA base-pair parameters\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "RNA base-pair step parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "RNA base-pair helical parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 1, heli_par, fp);

    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
  }

  /* Step and helical parameters along each strand */
  for (i = 1; i <= ds; i++) {
    for (j = 1; j <= nbpm1; j++) {
      refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
      helical_par(r1, o1, r2, o2, hpi, mfi, mfoi);
      if (i == 1) {
        pvt_dxdy(r1, o1, pvt1, hpi, mfi, mfoi);
        rna(r1, o1, pvt1, r2, o2, pvt1, step_par[j], mfi, mfoi);
      } else {
        pvt_dxdy(r1, o1, pvt2, hpi, mfi, mfoi);
        rna(r1, o1, pvt2, r2, o2, pvt2, step_par[j], mfi, mfoi);
      }
      for (k = 1; k <= 6; k++)
        heli_par[j][k] = hpi[k];
    }
    if (i == 1) {
      print_sep(fp, '*', 76);
      fprintf(fp, "Strand I base step parameters\n");
      print_par(bp_seq, num_bp, 3, 0, step_par, fp);
      print_sep(fp, '*', 76);
      fprintf(fp, "Strand I base helical parameters\n");
      print_par(bp_seq, num_bp, 3, 1, heli_par, fp);
    } else {
      print_sep(fp, '*', 76);
      fprintf(fp, "Strand II base step parameters\n");
      print_par(bp_seq, num_bp, 4, 0, step_par, fp);
      print_sep(fp, '*', 76);
      fprintf(fp, "Strand II base helical parameters\n");
      print_par(bp_seq, num_bp, 4, 1, heli_par, fp);
    }
  }

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mfi, 1, 3, 1, 3);
  free_dmatrix(step_par, 1, nbpm1, 1, 6);
  free_dmatrix(heli_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters for duplex based on all 7 methods */
void other_pars(long num_bp, char **bp_seq, long *bphlx, double **orien,
                double **org) {
  FILE *fp;

  fp = open_file(SEVEN_FILE, "w");

  out_cehs(num_bp, bp_seq, bphlx, orien, org, fp);
  out_compdna(num_bp, bp_seq, bphlx, orien, org, fp);
  out_curves(num_bp, bp_seq, bphlx, orien, org, fp);
  out_freehelix(num_bp, bp_seq, bphlx, orien, org, fp);
  out_ngeom(num_bp, bp_seq, bphlx, orien, org, fp);
  out_nuparm(num_bp, bp_seq, bphlx, orien, org, fp);
  out_rna(2L, num_bp, bp_seq, bphlx, orien, org, fp);

  close_file(fp);
}

/// @brief - app_fncs.cpp

/* global variables definition */
struct_Gvars Gvars;

static char *CMARKERS =
    "||x+0"; /* helix begin, middle, end, isolated bp, last */
static char *REBUILD_CIDS = "ABCDEFGHI";

static char *ptable[] = {
    " H", "HE", "LI", "BE", " B", " C", " N", " O", " F", "NE", "NA",
    "MG", "AL", "SI", " P", " S", "CL", "AR", " K", "CA", "SC", "TI",
    " V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS",
    "SE", "BR", "KR", "RB", "SR", " Y", "ZR", "NB", "MO", "TC", "RU",
    "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", " I", "XE", "CS",
    "BA", "LA", "HF", "TA", " W", "RE", "OS", "IR", "PT", "AU", "HG",
    "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "RF", "DB",
    "SG", "BH", "HS", "MT", "CE", "PR", "ND", "PM", "SM", "EU", "GD",
    "TB", "DY", "HO", "ER", "TM", "YB", "LU", "TH", "PA", " U", "NP",
    "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR"};

static void get_3dna_homedir(char *homedir) {
  char *temp;

  if ((temp = getenv("X3DNA")) != NULL) {
    strcpy(homedir, temp);
    check_slash(homedir); /* with ending slash */
    return;
  }

  fprintf(stderr, "Please set the X3DNA environment variable\n");
  fatal(
      "Run 'x3dna_setup', and visit http://forum.x3dna.org/ for more info.\n");
}

static void get_3dna_version(char *homedir, char *version) {
  char *p0, *line, str[BUF1K];
  FILE *fp;

  sprintf(str, "%sconfig/version", homedir);
  fp = open_file(str, "r");
  while ((p0 = my_getline(fp)) != NULL) {
    line = trim(p0);           /* keep the original address of p0 */
    if (!is_skip_line(line)) { /* just the first line */
      strcpy(version, line);
      free(p0);
      break;
    }
    free(p0);
  }
  close_file(fp);
}

int case_strcmp(const char *s1, const char *s2) {
  int i, c1, c2;

  for (i = 0; c1 = toupper((int)s1[i]), c2 = toupper((int)s2[i]), c1 == c2; i++)
    if (c1 == '\0')
      return 0;
  return c1 - c2;
}

int case_strncmp(const char *s1, const char *s2, long n) {
  int i, c1, c2;

  for (i = 0;
       (c1 = toupper((int)s1[i]), c2 = toupper((int)s2[i]), c1 == c2) && i < n;
       i++)
    if (c1 == '\0')
      return 0;

  return (i >= n) ? 0 : c1 - c2;
}

char *case_strstr(const char *haystack, const char *needle) {
  char *haystack_cp, *needle_cp, *p;

  haystack_cp = my_strdup(haystack);
  needle_cp = my_strdup(needle);

  lowerstr(haystack_cp);
  lowerstr(needle_cp);

  p = strstr(haystack_cp, needle_cp);

  if (p != NULL)
    p = (p - haystack_cp) + (char *)haystack;

  free_cvector(haystack_cp, 0, -1);
  free_cvector(needle_cp, 0, -1);

  return p;
}

char *case_strchr(const char *s, int c) {
  char *str, *p;

  str = my_strdup(s);
  lowerstr(str);

  p = strchr(str, tolower(c));

  if (p != NULL)
    p = (p - str) + (char *)s;

  free_cvector(str, 0, -1);

  return p;
}

long is_empty_string(const char *str) {
  return (str == NULL || strcmp(str, "") == 0);
}

long is_equal_string(const char *str1, const char *str2) {
  return (strcmp(str1, str2) == 0);
}

/* no processing of line[] here */
long is_skip_line(char *line) {
  if (strchr(SKIPS, *line)) /* line starting with # or empty */
    return true;
  else
    return false;
}

/* get a new name with base name from 'src' and extension 'ext' */
void bname_ext(char *src, char *ext, char *dst) {
  char str[BUF512];

  bname_noext(src, str);
  sprintf(dst, "%s.%s", str, ext);
}

double get_point2line_perp_distance(double *pnt, double *line_p1,
                                    double *line_p2) {
  long i;
  double d, p1_p2[4], p1_pnt[4];

  ddxyz(line_p1, line_p2, p1_p2);
  vec_norm(p1_p2);

  ddxyz(line_p1, pnt, p1_pnt);
  d = dot(p1_pnt, p1_p2);

  for (i = 1; i <= 3; i++)
    p1_pnt[i] -= d * p1_p2[i];

  return veclen(p1_pnt);
}

void get_tag_string_pair(char *prefix, char *tag, char *btag, char *etag) {
  sprintf(btag, "<%s%s", prefix, tag);
  sprintf(etag, "</%s%s", prefix, tag);
}

void get_xml_tag(FILE *fpxml, char *prefix, char *line, char *connector,
                 char *tag, char *tag_str) {
  char *p0, *p1, *ep, btag[BUF512], etag[BUF512], temp[BUFBIG];
  long isok = false;

  get_tag_string_pair(prefix, tag, btag, etag);

  if (!tag_match(prefix, line, tag))
    return;

  strcpy(tag_str, ""); /* initialize it */
  ep = strstr(line, etag);

  if (ep) { /* tag contained in one line */
    isok = true;
    *ep = '\0';
    strcpy(temp, strchr(line, '>') + 1);
    strcpy(tag_str, trim(temp));

  } else { /* in multiple lines */
    strcpy(tag_str, strchr(line, '>') + 1);
    while ((p0 = my_getline(fpxml)) != NULL) {
      p1 = trim(p0);
      ep = strstr(p1, etag);
      if (ep) { /* tag contained in one line */
        isok = true;
        *ep = '\0';
        strcat(tag_str, connector);
        strcat(tag_str, p1);
        free(p0);
        break;
      } else {
        strcat(tag_str, connector);
        strcat(tag_str, p1);
      }
      free(p0);
    }
  }

  if (!isok)
    fatal("tag <%s> has end match\n", tag);
}

void get_xml_tag_long(FILE *fpxml, char *prefix, char *line, char *tag,
                      long *lval) {
  char tag_str[BUFBIG];

  if (tag_match(prefix, line, tag)) {
    get_xml_tag(fpxml, prefix, line, " ", tag, tag_str);
    *lval = cvt2long(tag_str);
  }
}

void get_xml_tag_double(FILE *fpxml, char *prefix, char *line, char *tag,
                        double *dval) {
  char tag_str[BUFBIG];

  if (tag_match(prefix, line, tag)) {
    get_xml_tag(fpxml, prefix, line, " ", tag, tag_str);
    *dval = cvt2double(tag_str);
  }
}

long tag_match(char *prefix, char *line, char *tag) {
  char btag[BUF512];

  sprintf(btag, "<%s%s", prefix, tag);

  return str_pmatch(line, btag) &&
         !strstr(line, "/>"); /* not short-handed form */
}

long set_switch_default_true(char *option) {
  char str[BUF512];
  long switch_set = true;

  if (!strchr(option, '=')) /* just a switch */
    return switch_set;

  get_strvalue(option, str, false);
  lowerstr(str);

  if (str_pmatch(str, "of") ||
      str_pmatch(str, "0") || /* off | 0 | no | false */
      str_pmatch(str, "n") || str_pmatch(str, "f"))
    switch_set = false;

  return switch_set;
}

static void check_required_option(char *option, char *invalid_str, char *msg) {
  if (is_equal_string(option, invalid_str)) {
    fprintf(stderr, "missing required option: %s\n", msg);
    fatal("Please try \"%s -h\" for usage information\n", Gvars.PROGNAME);
  }
}

/* get the local reference frame for each base. only the ring atoms are
   included in least-squares fitting. similar to ref_frames() */
void base_frame(long num_residue, char *bseq, long **seidx, long *res_type,
                char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, char *BDIR, double **orien,
                double **org) {
  static char *RingAtom[] = {RA_LIST};
  char idmsg[BUF512], sidmsg[BUF512], spdb[BUF512];
  char *sChainID, **sAtomName, **sResName, **sMiscs;
  double **sxyz, **eRing_xyz, **sRing_xyz, **fitted_xyz, **R;
  long i, ib, ie, j, RingAtom_num, RR9 = 9;
  long exp_katom, std_katom, nmatch, snum, *sResSeq;

  eRing_xyz = dmatrix(1, RR9, 1, 3);
  sRing_xyz = dmatrix(1, RR9, 1, 3);

  sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
  sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
  sChainID = cvector(1, NUM_RESIDUE_ATOMS);
  sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
  sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
  sMiscs = cmatrix(1, NUM_RESIDUE_ATOMS, 0, NMISC);

  fitted_xyz = dmatrix(1, RR9, 1, 3);
  R = dmatrix(1, 3, 1, 3);

  for (i = 1; i <= num_residue; i++) {
    if (res_type[i] < 0)
      continue; /* non-bases */

    ib = seidx[i][1];
    ie = seidx[i][2];
    get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);

    RingAtom_num = (res_type[i] == 1) ? 9 : 6;
    set_std_base_pdb(BDIR, false, bseq[i], spdb);
    snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                    sMiscs, 1, "*");
    sprintf(sidmsg, "in standard base: %s", spdb);

    nmatch = 0;
    for (j = 0; j < RingAtom_num; j++) {
      exp_katom = find_1st_atom(RingAtom[j], AtomName, ib, ie, idmsg);
      std_katom = find_1st_atom(RingAtom[j], sAtomName, 1, snum, sidmsg);
      if (exp_katom && std_katom) {
        ++nmatch;
        cpxyz(xyz[exp_katom], eRing_xyz[nmatch]);
        cpxyz(sxyz[std_katom], sRing_xyz[nmatch]);
      }
    }

    (void)ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, org[i]);
    mst2orien(orien[i], 0, R);
  }

  free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq,
           sxyz, sMiscs);
  free_dmatrix(eRing_xyz, 1, RR9, 1, 3);
  free_dmatrix(sRing_xyz, 1, RR9, 1, 3);
  free_dmatrix(fitted_xyz, 1, RR9, 1, 3);
  free_dmatrix(R, 1, 3, 1, 3);
}

static void align_block_xyz(double *orien_i, double *org_i, double *xyz_j,
                            double *new_xyz) {
  long i, j;
  double temp[4];

  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++)
      temp[j] = orien_i[i + (j - 1) * 3];
    new_xyz[i] = dot(xyz_j, temp) + org_i[i];
  }
}

static void reset_alist_symbol_idx(char *alist, long *aidx) {
  char c, asym[3], atoms_list[NELE][3];
  long i, k, nchar, natom;

  nchar = (long)strlen(alist);
  if (nchar % 2 != 0)
    fatal("each atomic symbol must be two-char long <%s>\n", alist);
  if (nchar == 0)
    fatal("empty atom list not allowed\n");

  upperstr(alist); /* to all upper case */
  for (i = 0; i < nchar; i++) {
    c = alist[i];
    if ((c != '.') && !isupper((int)c))
      fatal("invalid char '%c': must be [.A-Z]\n", c);
    if (c == '.')
      alist[i] = ' ';
  }

  natom = nchar / 2;
  aidx[0] = natom; /* number of atoms in the list */

  atom_info(1, atoms_list, NULL, NULL);

  for (i = 1; i <= natom; i++) {
    k = (i - 1) * 2; /* 0, 2, 4 ... index position for next atom symbol */
    strncpy(asym, alist + k, 2);
    asym[2] = '\0';

    if (!num_strmatch(asym, Gvars.ATOM_NAMES, 0, Gvars.NUM_ELE))
      fatal("invalid atom symbol <%s>\n", asym);

    aidx[i] = asym_idx(asym, atoms_list, 0);
  }
}

static void check_space_in_altlist(char *alt_list) {
  char str[BUF1K];

  if (strchr(alt_list, ' ') == NULL) { /* not included yet */
    sprintf(str, " %s", alt_list);
    strcpy(alt_list, str);
  }
}

static void check_if_valid_xml_parfile(void) {
  static char *pars[] = {
      "min_base_hb", "hb_lower",   "hb_dist1",        "hb_dist2",
      "hb_atoms",    "alt_list",   "max_dorg",        "min_dorg",
      "max_dv",      "min_dv",     "max_plane_angle", "min_plane_angle",
      "max_dNN",     "min_dNN",    "helix_break",     "std_curved",
      "water_dist",  "water_dlow", "water_atoms",     "o3p_dist"};
  char *prefix = "", *dft = "xxxxxx", str[BUF512];
  char *p0, *line, fileLoc[BUF512];
  long i, num_pars = sizeof pars / sizeof pars[0];
  FILE *fp;

  get_BDIR(fileLoc, PAR_FILE);
  strcat(fileLoc, PAR_FILE);

  fp = open_file(fileLoc, "r");

  while ((p0 = my_getline(fp)) != NULL) {
    line = trim(p0); /* keep the original value of p0 */
    if (is_skip_line(line)) {
      free(p0);
      continue;
    }
    for (i = 0; i < num_pars; i++) {
      strcpy(str, dft);
      get_xml_tag(fp, prefix, line, " ", pars[i], str);
      if (!is_equal_string(str, dft)) { /* with at least one valid xml tag */
        free(p0);
        close_file(fp);
        return;
      }
    }
    free(p0);
  }

  close_file(fp);

  fprintf(stderr,
          "File <%s> does NOT contain any valid xml tags for parameters\n"
          "\tmaybe it is in the previous v1.5 version?\n"
          "\tso system default setting is used here...\n",
          fileLoc);
}

static void mask_off_options(miscPars *misc_pars) {
  misc_pars->min_base_hb = 1;
  misc_pars->hb_lower = 0.0;
  misc_pars->min_dorg = 0.0;
  misc_pars->min_dv = 0.0;
  misc_pars->min_plane_angle = 0.0;
  misc_pars->max_dNN = XBIG;
  misc_pars->water_dlow = 0.0;
}

void set_default_misc_pars(miscPars *misc_pars) {
  misc_pars->min_base_hb = 1; /* at least one H-bond between base atoms */
  /* 1jgq/1jgo/1jgp/1gix has N2 * N4  1.77 for pair C1_G-**+-C_C2 */
  misc_pars->hb_lower = 1.8;
  misc_pars->hb_dist1 = 4.0;
  misc_pars->hb_dist2 = 0.0;
  strcpy(misc_pars->hb_atoms, ".O.N");
  strcpy(misc_pars->alt_list, "A1");

  misc_pars->max_dorg = 15.0;
  misc_pars->min_dorg = 0.0;
  misc_pars->max_dv = 2.5;
  misc_pars->min_dv = 0.0;
  misc_pars->max_plane_angle = 65.0;
  misc_pars->min_plane_angle = 0.0;
  misc_pars->max_dNN = XBIG;
  misc_pars->min_dNN = 4.5;

  misc_pars->helix_break = 7.5;
  misc_pars->std_curved = 0.6;

  misc_pars->water_dist = 3.2;
  misc_pars->water_dlow = 0.0;
  strcpy(misc_pars->water_atoms, ".O.N");

  misc_pars->o3p_dist = 4.5;
}

static void output_misc_pars(FILE *fp, miscPars *misc_pars, char sep) {
  long k = 76;

  print_sep(fp, sep, k);

  fprintf(fp, "<min_base_hb>%ld</min_base_hb>\n", misc_pars->min_base_hb);
  fprintf(fp, "<hb_lower>%g</hb_lower>\n", misc_pars->hb_lower);
  fprintf(fp, "<hb_dist1>%g</hb_dist1>\n", misc_pars->hb_dist1);
  fprintf(fp, "<hb_dist2>%g</hb_dist2>\n", misc_pars->hb_dist2);
  fprintf(fp, "<hb_atoms>%s</hb_atoms>\n", misc_pars->hb_atoms);
  fprintf(fp, "<alt_list>%s</alt_list>\n", misc_pars->alt_list);

  fprintf(fp, "<max_dorg>%g</max_dorg>\n", misc_pars->max_dorg);
  fprintf(fp, "<min_dorg>%g</min_dorg>\n", misc_pars->min_dorg);
  fprintf(fp, "<max_dv>%g</max_dv>\n", misc_pars->max_dv);
  fprintf(fp, "<min_dv>%g</min_dv>\n", misc_pars->min_dv);
  fprintf(fp, "<max_plane_angle>%g</max_plane_angle>\n",
          misc_pars->max_plane_angle);
  fprintf(fp, "<min_plane_angle>%g</min_plane_angle>\n",
          misc_pars->min_plane_angle);
  fprintf(fp, "<max_dNN>%g</max_dNN>\n", misc_pars->max_dNN);
  fprintf(fp, "<min_dNN>%g</min_dNN>\n", misc_pars->min_dNN);

  fprintf(fp, "<helix_break>%g</helix_break>\n", misc_pars->helix_break);
  fprintf(fp, "<std_curved>%g</std_curved>\n", misc_pars->std_curved);

  fprintf(fp, "<water_dist>%g</water_dist>\n", misc_pars->water_dist);
  fprintf(fp, "<water_dlow>%g</water_dlow>\n", misc_pars->water_dlow);
  fprintf(fp, "<water_atoms>%s</water_atoms>\n", misc_pars->water_atoms);

  fprintf(fp, "<o3p_dist>%g</o3p_dist>\n", misc_pars->o3p_dist);

  print_sep(fp, sep, k);
  fprintf(fp, "\n");
}

static void check_misc_pars(miscPars *misc_pars) {
  if (misc_pars->min_base_hb < 0 || misc_pars->min_base_hb > 3) {
    fprintf(stderr, "min_base_hb: <%ld> -- reset to the default (1)\n",
            misc_pars->min_base_hb);
    misc_pars->min_base_hb = 1;
  }

  if (misc_pars->hb_lower < 0.0)
    fatal("hb_lower (%g) < 0.0\n", misc_pars->hb_lower);

  if (misc_pars->hb_dist1 < 0.0)
    fatal("hb_dist1 (%g) < 0.0\n", misc_pars->hb_dist1);

  if (misc_pars->hb_dist2 < 0.0)
    fatal("hb_dist2 (%g) < 0.0\n", misc_pars->hb_dist2);

  reset_alist_symbol_idx(misc_pars->hb_atoms, misc_pars->hb_idx);
  check_space_in_altlist(misc_pars->alt_list);

  if (misc_pars->max_dorg < 0.0)
    fatal("max_dorg (%g) < 0.0\n", misc_pars->max_dorg);

  if (misc_pars->min_dorg < 0.0)
    fatal("min_dorg (%g) < 0.0\n", misc_pars->min_dorg);
  if (misc_pars->min_dorg > misc_pars->max_dorg)
    fatal("min_dorg (%g) > max_dorg (%g)\n", misc_pars->min_dorg,
          misc_pars->max_dorg);

  if (misc_pars->max_dv < 0.0)
    fatal("max_dv (%g) < 0.0\n", misc_pars->max_dv);

  if (misc_pars->min_dv < 0.0)
    fatal("min_dv (%g) < 0.0\n", misc_pars->min_dv);
  if (misc_pars->min_dv > misc_pars->max_dv)
    fatal("min_dv (%g) > max_dv (%g)\n", misc_pars->min_dv, misc_pars->max_dv);

  if ((misc_pars->max_plane_angle < 0.0) || (misc_pars->max_plane_angle > 90.0))
    fatal("max_plane_angle (%g) must be in range [0, 90]\n",
          misc_pars->max_plane_angle);

  if ((misc_pars->min_plane_angle < 0.0) || (misc_pars->min_plane_angle > 90.0))
    fatal("min_plane_angle (%g) must be in range [0, 90]\n",
          misc_pars->min_plane_angle);
  if (misc_pars->min_plane_angle > misc_pars->max_plane_angle)
    fatal("min_plane_angle (%g) > max_plane_angle (%g)\n",
          misc_pars->min_plane_angle, misc_pars->max_plane_angle);

  if (misc_pars->max_dNN < 0.0)
    fatal("max_dNN (%g) < 0.0\n", misc_pars->max_dNN);

  if (misc_pars->min_dNN < 0.0)
    fatal("min_dNN (%g) < 0.0\n", misc_pars->min_dNN);
  if (misc_pars->min_dNN > misc_pars->max_dNN)
    fatal("min_dNN (%g) > max_dNN (%g)\n", misc_pars->min_dNN,
          misc_pars->max_dNN);

  if (misc_pars->helix_break < 0.0)
    fatal("helix_break (%g) < 0.0\n", misc_pars->helix_break);

  if (misc_pars->std_curved < 0.0)
    fatal("std_curved (%g) < 0.0\n", misc_pars->std_curved);

  if (misc_pars->water_dist < 0.0)
    fatal("water_dist (%g) < 0.0\n", misc_pars->water_dist);

  if (misc_pars->water_dlow < 0.0)
    fatal("water_dlow (%g) < 0.0\n", misc_pars->water_dlow);
  if (misc_pars->water_dlow > misc_pars->water_dist)
    fatal("water_dlow (%g) > water_dist (%g)\n", misc_pars->water_dlow,
          misc_pars->water_dist);

  reset_alist_symbol_idx(misc_pars->water_atoms, misc_pars->water_idx);

  if (misc_pars->o3p_dist < 0.0)
    fatal("o3p_dist (%g) < 0.0\n", misc_pars->o3p_dist);
}

static void get_3dna_pars(miscPars *misc_pars) {
  char *p0, *line, *prefix = "", fileLoc[BUF512];
  FILE *fp, *fp2 = NULL;

  check_if_valid_xml_parfile();
  set_default_misc_pars(misc_pars);

  if (Gvars.DEBUG >= DEBUG_LEVEL) {
    fp2 = open_file("chk_misc_3dna.par", "w");
    output_misc_pars(fp2, misc_pars, '0');
  }

  get_BDIR(fileLoc, PAR_FILE);
  strcat(fileLoc, PAR_FILE);
  if (Gvars.VERBOSE)
    fprintf(stderr, " ...... reading file: %s ...... \n", PAR_FILE);

  fp = open_file(fileLoc, "r");

  while ((p0 = my_getline(fp)) != NULL) {
    line = trim(p0); /* keep the original value of p0 */
    if (is_skip_line(line)) {
      free(p0);
      continue;
    }

    get_xml_tag_long(fp, prefix, line, "min_base_hb", &misc_pars->min_base_hb);
    get_xml_tag_double(fp, prefix, line, "hb_lower", &misc_pars->hb_lower);
    get_xml_tag_double(fp, prefix, line, "hb_dist1", &misc_pars->hb_dist1);
    get_xml_tag_double(fp, prefix, line, "hb_dist2", &misc_pars->hb_dist2);
    get_xml_tag(fp, prefix, line, " ", "hb_atoms", misc_pars->hb_atoms);
    get_xml_tag(fp, prefix, line, " ", "alt_list", misc_pars->alt_list);

    get_xml_tag_double(fp, prefix, line, "max_dorg", &misc_pars->max_dorg);
    get_xml_tag_double(fp, prefix, line, "min_dorg", &misc_pars->min_dorg);
    get_xml_tag_double(fp, prefix, line, "max_dv", &misc_pars->max_dv);
    get_xml_tag_double(fp, prefix, line, "min_dv", &misc_pars->min_dv);
    get_xml_tag_double(fp, prefix, line, "max_plane_angle",
                       &misc_pars->max_plane_angle);
    get_xml_tag_double(fp, prefix, line, "min_plane_angle",
                       &misc_pars->min_plane_angle);
    get_xml_tag_double(fp, prefix, line, "max_dNN", &misc_pars->max_dNN);
    get_xml_tag_double(fp, prefix, line, "min_dNN", &misc_pars->min_dNN);

    get_xml_tag_double(fp, prefix, line, "helix_break",
                       &misc_pars->helix_break);
    get_xml_tag_double(fp, prefix, line, "std_curved", &misc_pars->std_curved);

    get_xml_tag_double(fp, prefix, line, "water_dist", &misc_pars->water_dist);
    get_xml_tag_double(fp, prefix, line, "water_dlow", &misc_pars->water_dlow);
    get_xml_tag(fp, prefix, line, " ", "water_atoms", misc_pars->water_atoms);

    get_xml_tag_double(fp, prefix, line, "o3p_dist", &misc_pars->o3p_dist);

    free(p0);
  }

  close_file(fp);

  if (Gvars.DEBUG >= DEBUG_LEVEL)
    output_misc_pars(fp2, misc_pars, '1');

  check_misc_pars(misc_pars);

  if (Gvars.DEBUG >= DEBUG_LEVEL) {
    output_misc_pars(fp2, misc_pars, '2');
    close_file(fp2);
  }

  if (Gvars.DEBUG < DEBUG_LEVEL)
    mask_off_options(misc_pars);
}

static long overwrite_misc_pars(char *option) {
  miscPars *misc_pars = &Gvars.misc_pars;

  if (case_str_pmatch(option, "-min_base_hb")) {
    misc_pars->min_base_hb = get_lvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-hb_lower")) {
    misc_pars->hb_lower = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-hb_dist1")) {
    misc_pars->hb_dist1 = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-hb_dist2")) {
    misc_pars->hb_dist2 = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-hb_atoms")) {
    get_strvalue(option, misc_pars->hb_atoms, false);
    reset_alist_symbol_idx(misc_pars->hb_atoms, misc_pars->hb_idx);
    return true;

  } else if (case_str_pmatch(option, "-alt_list")) {
    get_strvalue(option, misc_pars->alt_list, false);
    check_space_in_altlist(misc_pars->alt_list);
    return true;

  } else if (case_str_pmatch(option, "-max_dorg")) {
    misc_pars->max_dorg = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-min_dorg")) {
    misc_pars->min_dorg = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-max_dv")) {
    misc_pars->max_dv = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-min_dv")) {
    misc_pars->min_dv = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-max_plane_angle")) {
    misc_pars->max_plane_angle = get_dvalue(option, 0, 90);
    return true;

  } else if (case_str_pmatch(option, "-min_plane_angle")) {
    misc_pars->min_plane_angle = get_dvalue(option, 0, 90);
    return true;

  } else if (case_str_pmatch(option, "-max_dNN")) {
    misc_pars->max_dNN = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-min_dNN")) {
    misc_pars->min_dNN = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-helix_break")) {
    misc_pars->helix_break = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-std_curved")) {
    misc_pars->std_curved = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-water_dist")) {
    misc_pars->water_dist = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-water_dlow")) {
    misc_pars->water_dlow = get_dvalue(option, 0, 999);
    return true;

  } else if (case_str_pmatch(option, "-water_atoms")) {
    get_strvalue(option, misc_pars->water_atoms, false);
    reset_alist_symbol_idx(misc_pars->water_atoms, misc_pars->water_idx);
    return true;

  } else if (case_str_pmatch(option, "-o3p_dist")) {
    misc_pars->o3p_dist = get_dvalue(option, 0, 999);
    return true;
  }

  return false;
}

/* base-pairing type information as in check_pair():
 * +1 WC geometry; +2: WC pair; -1: other cases */
long read_PairInfo(char *inpfile, long **pair_info) {
  char *p0, *p;
  long i, num_bp;
  FILE *fp;

  fp = open_file(inpfile, "r");

  for (i = 1; i <= 5; i++) { /* skip top 5 lines */
    p0 = my_getline(fp);
    if (i == 4 && sscanf(p0, "%ld", &num_bp) != 1)
      fatal("\tcan't extract number of base-pairs <%s>\n", p0);
    free(p0);
  }

  for (i = 1; i <= num_bp; i++) {
    p0 = my_getline(fp);
    if (sscanf(p0, "%ld %ld", &pair_info[i][1], &pair_info[i][2]) != 2)
      fatal("\tcan't extract pairing number %ld <%s>\n", i, p0);
    pair_info[i][0] = -1;
    p = strchr(p0, ']') + 3; /* ]G-----C[ */
    if (*(p + 1) == '-')     /* WC geometry */
      pair_info[i][0] = 1;
    if (*p == '-') /* WC pair */
      pair_info[i][0] = 2;
    free(p0);
  }

  close_file(fp);

  return num_bp;
}

/* use O3' and P distance to check if residues i & j are directly connected */
long is_linked(long i, long j, double **o3_p) {
  double d;

  d = distance_ab(o3_p, i, j, 4, 8); /* O3'[i]-P[j] */
  if (dval_in_range(d, 0.0, O3P_UPPER))
    return 1L; /* 5'--->3' linkage */

  d = distance_ab(o3_p, j, i, 4, 8); /* O3'[j]-P[i] */
  if (dval_in_range(d, 0.0, O3P_UPPER))
    return -1L; /* 3'--->5' linkage */

  return 0L; /* no direct linkage */
}

/* calculate distance between ia & ib: ipa & ipb mark if ia/ib exist */
double distance_ab(double **o3_p, long ia, long ib, long ipa, long ipb) {
  return (o3_p[ia][ipa] > 0.0 && o3_p[ib][ipb] > 0.0)
             ? p1p2_dist(o3_p[ia] + ipa - 4, o3_p[ib] + ipb - 4)
             : -1.0;
}

void reverse_y_z_columns(double **R) {
  long i;

  for (i = 1; i <= 3; i++) {
    R[i][2] = -R[i][2]; /* reverse y-axis */
    R[i][3] = -R[i][3]; /* reverse z-axis */
  }
}

/* this function must be run before calling atom_idx() */
void set_my_globals(char *pgname) {

  Gvars.VERBOSE = false;
  Gvars.CHAIN_CASE = true; /* case-sensitive */
  Gvars.ALL_MODEL = false; /* just one model: ENDMDL/END */
  Gvars.DEBUG = false;     /* for general distribution */
  Gvars.PROGNAME = pgname;
  Gvars.ATTACH_RESIDUE = true;    /* add connected HETATM metals or residues */
  Gvars.THREE_LETTER_NTS = false; /* ADE/CYT/GUA/THY/URA */
  Gvars.PDBV3 = true;
  Gvars.ORIGINAL_COORDINATE = false;
  Gvars.OCCUPANCY = false;
  Gvars.HEADER = true;

  Gvars.mmcif = false;
  Gvars.NT_CUTOFF = 0.2618; /* 2o8b C30 on F has 0.24 */

  get_3dna_homedir(Gvars.X3DNA_HOMEDIR);
  get_3dna_version(Gvars.X3DNA_HOMEDIR, Gvars.X3DNA_VER);
  strcpy(Gvars.CHAIN_MARKERS, CMARKERS);
  strcpy(Gvars.REBUILD_CHAIN_IDS, REBUILD_CIDS);

  Gvars.ATOM_NAMES = ptable;
  Gvars.NUM_ELE = sizeof ptable / sizeof ptable[0] - 1; /* minus - 1 */

  Gvars.ATOMLIST = cmatrix(1, BUFBIG, 0, 6);
  get_atomlist(Gvars.ATOMLIST, &Gvars.NUM_SATOM);

  Gvars.BASELIST = cmatrix(1, BUFBIG, 0, 4);
  get_baselist(Gvars.BASELIST, &Gvars.NUM_SBASE);

  Gvars.AtomName0 = NULL;
  Gvars.ResName0 = NULL;
  Gvars.Name0 = false;

  Gvars.label_RC8_YC6 = false;

  get_3dna_pars(&Gvars.misc_pars);
}

void clear_my_globals(void) {
  free_cmatrix(Gvars.ATOMLIST, 1, -1, 0, -1);
  free_cmatrix(Gvars.BASELIST, 1, -1, 0, -1);
}

static void set_chain_markers(char *option) {
  long k;

  get_strvalue(option, Gvars.CHAIN_MARKERS, false);
  k = strlen(Gvars.CHAIN_MARKERS);
  if (k < 4) {
    fprintf(stderr, "too short [< 4-char] input for option: %s\n", option);
    fprintf(stderr, "\treset to default: -chain_markers='%s'\n", CMARKERS);
    strcpy(Gvars.CHAIN_MARKERS, CMARKERS);
  } else if (k == 4) {
    if (Gvars.CHAIN_MARKERS[0] == CMARKERS[0])
      Gvars.CHAIN_MARKERS[4] = CMARKERS[4];
    else
      Gvars.CHAIN_MARKERS[4] = '1';
  }
}

static void set_rebuild_cids(char *option) {
  char cids[BUF512];
  long k;

  strcpy(cids, "");
  get_strvalue(option, cids, false);
  k = strlen(cids);

  if (k == 0)
    strcpy(cids, REBUILD_CIDS);
  else if (k == 1)
    strcat(cids, "B");

  cids[2] = '\0'; /* keep only 2 chars */
  fprintf(stderr, "[i] -- chain ids for rebuild: '%s'\n", cids);

  strcpy(Gvars.REBUILD_CHAIN_IDS, cids);
}

long check_global_options(char *option) {
  if (lux_ncmatch(option, "^--?debug")) {
    Gvars.DEBUG = get_lvalue(option, 0, BUF512);
    return true;

  } else if (lux_ncmatch(option, "^--?chain_case")) {
    Gvars.CHAIN_CASE = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?chain_markers")) {
    set_chain_markers(option);
    return true;

  } else if (lux_ncmatch(option, "^--?rebuild-cids")) {
    set_rebuild_cids(option);
    return true;

  } else if (lux_ncmatch(option, "^--?(all_model|models|symm)")) {
    Gvars.ALL_MODEL = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?verbose")) {
    Gvars.VERBOSE = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?attach")) {
    Gvars.ATTACH_RESIDUE = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?three")) {
    Gvars.THREE_LETTER_NTS = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?pdbv3")) {
    Gvars.PDBV3 = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?(ori|raw)")) {
    Gvars.ORIGINAL_COORDINATE = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?occ")) {
    Gvars.OCCUPANCY = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?hea")) {
    Gvars.HEADER = set_switch_default_true(option);
    if (Gvars.HEADER && lux_ncmatch(option, "^--?(whole|all"))
      Gvars.HEADER = BUF32;
    return true;

  } else if (lux_ncmatch(option, "^--?(mm)?cif")) {
    Gvars.mmcif = set_switch_default_true(option);
    return true;

  } else if (lux_ncmatch(option, "^--?ntc")) {
    Gvars.NT_CUTOFF = get_dvalue(option, 0.05, 1.0);
    return true;

  } else if (lux_ncmatch(option, "^--?label[-_]?C")) {
    Gvars.label_RC8_YC6 = set_switch_default_true(option);
    return true;

  } else if (overwrite_misc_pars(option))
    return true;

  else
    return false;
}

void residue_strid(char chain_id, long res_seq, char *misc, char *rname,
                   char *idmsg) {
  char modelNum[10], iCode = misc[2], newName[BUF512];
  long model_lval, nlen = 4;

  strncpy(modelNum, misc + 30, nlen);
  modelNum[nlen] = '\0';
  model_lval = cvt2long(modelNum); /* model number of the residue */

  if (iCode == ' ')
    iCode = '.';
  convert_resNameSpace(rname, '_', newName); /* '\0' to eliminate space */

  if (model_lval == 0)
    sprintf(idmsg, "%c%ld%c%s", chain_id, res_seq, iCode, newName);
  else
    sprintf(idmsg, "%4s%c%ld%c%s", modelNum, chain_id, res_seq, iCode, newName);
}

void convert_resNameSpace(char *resName, char replacement, char *newName) {
  long i, j = 0, k = strlen(resName);

  if (replacement != '\0') {
    for (i = 0; i < k; i++)
      newName[i] = (resName[i] == ' ') ? replacement : resName[i];
    newName[i] = '\0';
  } else { /* eliminate ' ' */
    for (i = 0; i < k; i++)
      if (resName[i] != ' ')
        newName[j++] = resName[i];
    newName[j] = '\0';
  }
}

void snap_atype(char **AtomName, long num_residue, long **seidx, long *res_type,
                long **atom_cidx) {
  char *BASE_ATOMS[] = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ",
                        " N7 ", " C8 ", " N9 ", " N6 ", " O2 ", " O4 ",
                        " C5M", " O6 ", " N2 ", " N4 "};
  char *PRBB_ATOMS[] = {" N  ", " C  ", " O  "};
  long num_dna = sizeof BASE_ATOMS / sizeof BASE_ATOMS[0] - 1;
  long num_pro = sizeof PRBB_ATOMS / sizeof PRBB_ATOMS[0] - 1;
  long i, ib, ie, j;

  for (i = 1; i <= num_residue; i++) {
    if (res_type[i] == -2)
      continue;
    ib = seidx[i][1];
    ie = seidx[i][2];
    for (j = ib; j <= ie; j++) {
      if (atom_cidx[1][j] == 3) /* skip H atoms */
        continue;
      if (res_type[i] == -1) { /* amino-acid */
        if (num_strmatch(AtomName[j], PRBB_ATOMS, 0, num_pro))
          atom_cidx[2][j] = -WITH_BKBN;
        else /* side chain */
          atom_cidx[2][j] = -WITH_BASE;
      } else { /* nucleotide */
        if (num_strmatch(AtomName[j], BASE_ATOMS, 0, num_dna))
          atom_cidx[2][j] = WITH_BASE;
        else
          atom_cidx[2][j] = WITH_BKBN;
      }
    }
  }
}

void cleanup_files(long renew, long cleanup) {
  static char *SAA[] = {AA_LIST};
  static char *RY[] = {"AT", "GC"};
  static char *DNA = DNA_BASE;
  static char *ALC[] = {"atom_lkg", "ablk_lkg", "bblk_lkg", "pblk_lkg",
                        "pall_lkg", "snap_all", "snap_blk"};
  static char *MISC[] = {"snap_bp.dat", "snap_nts.pdb", "snap_options"};
  char filename[BUF512], suffix[BUF512];
  long i, j, k, num_saa, num_ry, num_dna, num_alc, num_misc;

  if (!renew && !cleanup)
    return;

  num_saa = sizeof SAA / sizeof SAA[0] - 1;
  num_ry = sizeof RY / sizeof RY[0] - 1;
  num_dna = strlen(DNA) - 1;
  num_alc = sizeof ALC / sizeof ALC[0] - 1;
  num_misc = sizeof MISC / sizeof MISC[0] - 1;

  for (i = 0; i <= num_dna; i++) {
    for (j = 0; j <= num_saa; j++) {
      for (k = -1; k <= 3; k++) {
        if (k < 0) /* for backward compatibility */
          strcpy(suffix, "");
        else
          sprintf(suffix, "_%ld", k);
        sprintf(filename, "%c-%s%s.par", DNA[i], SAA[j], suffix);
        remove_file(filename);
        sprintf(filename, "%c-%s%s.pdb", DNA[i], SAA[j], suffix);
        remove_file(filename);
      }
    }
  }

  for (i = 0; i <= num_ry; i++) {
    for (j = 0; j <= num_saa; j++) {
      for (k = -1; k <= 3; k++) {
        if (k < 0) /* for backward compatibility */
          strcpy(suffix, "");
        else
          sprintf(suffix, "_%ld", k);
        sprintf(filename, "%s-%s%s.par", RY[i], SAA[j], suffix);
        remove_file(filename);
        sprintf(filename, "%s-%s%s.pdb", RY[i], SAA[j], suffix);
        remove_file(filename);
      }
    }
  }

  for (i = 0; i <= num_alc; i++) {
    sprintf(filename, "%s.alc", ALC[i]);
    remove_file(filename);
  }

  for (i = 0; i <= num_misc; i++)
    remove_file(MISC[i]);

  if (cleanup)
    fatal("done with cleaning up snap/poco files.\n");
}

void set_std_base_pdb(char *bdir, long irna, char bname, char *spdb) {
  char mb[BUF32], str[BUF1K];

  UNUSED_PARAMETER(bdir); /* for compatibility */

  strcpy(mb, irna ? "r" : "");

  /* current directory first */
  if (isupper((int)bname))
    sprintf(spdb, "%sAtomic_%c.pdb", mb, bname);
  else /* for modified bases */
    sprintf(spdb, "%sAtomic.%c.pdb", mb, bname);
  if (exist_file(spdb))
    return;

  /* the system directory */
  sprintf(str, "%sconfig/%s", Gvars.X3DNA_HOMEDIR, spdb);
  strcpy(spdb, str);
}

void parcat(char *str, double par, char *format, char *bstr) {
  char temp[BUF512];

  if (par > EMPTY_CRITERION) {
    sprintf(temp, format, par);
    strcat(str, temp);
  } else
    strcat(str, bstr);
}

/* endofline: check for and consume \r, \n, \r\n, or EOF */
static int endofline(FILE *fp, int c) {
  int eol;

  eol = (c == '\r' || c == '\n');
  if (c == '\r') {
    c = getc(fp);
    if (c != '\n' && c != EOF)
      ungetc(c, fp); /* read too far; put c back */
  }
  return eol;
}

static char *enlarge_cline(long *maxline, char *line)
/* for a 0-index char-array */
{
  char *newline;

  *maxline *= SFACTOR;

  if ((newline = (char *)realloc(line, (*maxline) * sizeof(char))) == NULL)
    fatal("realloc failure in enlarge_cline()\n");

  return newline;
}

char *my_getline(FILE *fp)
/* read line of arbitrary length from fp */
{
  int c;
  long i, maxline = BUF512;
  char *line = NULL;

  line = (char *)malloc(maxline * sizeof(char));
  if (line == NULL)
    fatal("out of memory in my_getline()\n");

  for (i = 0; (c = getc(fp)) != EOF && !endofline(fp, c); i++) {
    if (i >= maxline - 1)
      line = enlarge_cline(&maxline, line);
    line[i] = c;
  }
  line[i] = '\0';

  if (c == EOF && i == 0) {
    free(line); /* to avoid memory leak here */
    line = NULL;
  }

  return line;
}

/* trim leading and trailing white spaces */
char *trim(char *a) {
  int c;

  while (isspace(c = *a))
    a++;
  for (c = strlen(a) - 1; c >= 0; c--)
    if (!isspace((int)a[c]))
      break;
  a[c + 1] = '\0';

  return a;
}

/* itemize a string into its components */
long itemize(char *str, char *item[], long itemsize) {
  static char *sep_chars = " \t\r\n";
  char *p;
  long nitem = 0;

  for (p = strtok(str, sep_chars); p != NULL; p = strtok(NULL, sep_chars))
    if (nitem > itemsize) {
      fprintf(stderr, "ignoring items following %ldth\n", nitem);
      return itemsize;
    } else
      item[nitem++] = p;
  if (!nitem)
    fatal("empty input line not allowed\n");
  return nitem - 1; /* base-pair does not count */
}

/* itemize a string into tokens, and return number of token fields: 1-index */
long item_list(char *str, char *item[], long itemsize, char *sep_chars) {
  char *p;
  long nitem = 0;

  for (p = strtok(str, sep_chars); p != NULL; p = strtok(NULL, sep_chars)) {
    item[++nitem] = trim(p); /* get rid of leading and trailing spaces */
    if (nitem >= itemsize)
      return itemsize;
  }

  return nitem;
}

/* right-left base in a pair */
void refs_right_left(long bnum, double **orien, double **org, double **r1,
                     double *o1, double **r2, double *o2) {
  long ioffset3, ioffset9;

  ioffset3 = (bnum - 1) * 3;
  ioffset9 = (bnum - 1) * 9;

  cpxyz(org[2] + ioffset3, o1);
  cpxyz(org[1] + ioffset3, o2);
  orien2mst(orien[2], ioffset9, r1);
  orien2mst(orien[1], ioffset9, r2);
}

/* base (-pair) step */
void refs_i_j(long b1, long b2, double *bp_orien, double *bp_org, double **r1,
              double *o1, double **r2, double *o2) {
  ref_frame_i(b1, bp_orien, bp_org, r1, o1);
  ref_frame_i(b2, bp_orien, bp_org, r2, o2);
}

/* get the reference frame for "bnum" */
void ref_frame_i(long bnum, double *bp_orien, double *bp_org, double **r,
                 double *o) {
  long ioffset3, ioffset9;

  ioffset3 = (bnum - 1) * 3;
  ioffset9 = (bnum - 1) * 9;

  cpxyz(bp_org + ioffset3, o);
  orien2mst(bp_orien, ioffset9, r);
}

/* set orientation vector with ioffset from a 3-by-3 column-wise matrix */
void mst2orien(double *orien_vec, long ioffset, double **mst) {
  long i, ik, j;

  for (i = 1; i <= 3; i++) {
    ik = ioffset + (i - 1) * 3;
    for (j = 1; j <= 3; j++)
      orien_vec[ik + j] = mst[j][i];
  }
}

/* get 3-by-3 column-wise matrix from orientation vector with ioffset */
void orien2mst(double *orien_vec, long ioffset, double **mst) {
  long i, ik, j;

  for (i = 1; i <= 3; i++) {
    ik = ioffset + (i - 1) * 3;
    for (j = 1; j <= 3; j++)
      mst[j][i] = orien_vec[ik + j];
  }
}

/* given 3 vectors (x, y & z), get matrix (column-wise) */
void x_y_z_2_mtx(double *x, double *y, double *z, double **mtx) {
  long i;

  for (i = 1; i <= 3; i++) {
    mtx[i][1] = x[i];
    mtx[i][2] = y[i];
    mtx[i][3] = z[i];
  }
}

/* given a matrix (column-wise), get 3 vectors (x, y & z) */
void mtx_2_x_y_z(double **mtx, double *x, double *y, double *z) {
  long i;

  for (i = 1; i <= 3; i++) {
    x[i] = mtx[i][1];
    y[i] = mtx[i][2];
    z[i] = mtx[i][3];
  }
}

/* Averaging method 1 (CEHS way): always used for a single base-pair */
void cehs_average(long inum_base, long *ivec, double **orien, double **org,
                  double **mst, double *morg) {
  double pars[7], orgn[4], **bi, **mstn;
  long i, ik;

  mstn = dmatrix(1, 3, 1, 3);
  bi = dmatrix(1, 3, 1, 3);

  ik = ivec[1]; /* middle frame initialized to 1st base */
  cpxyz(org[ik], morg);
  orien2mst(orien[ik], 0, mst);
  for (i = 1; i <= inum_base; i++) { /* for each residue */
    ik = ivec[i];
    orien2mst(orien[ik], 0, bi);
    if (dot(&orien[ik][6], &orien[ivec[1]][6]) < 0.0)
      reverse_y_z_columns(bi);
    bpstep_par(bi, org[ik], mst, morg, pars, mstn, orgn);
    copy_dmatrix(mstn, 3, 3, mst);
    cpxyz(orgn, morg);
  }

  free_dmatrix(mstn, 1, 3, 1, 3);
  free_dmatrix(bi, 1, 3, 1, 3);
}

/* Write a base-pair or a multiplet with reference to its middle-frame:
 * There are 2 methods to calculate the middle-frame with > 2 bases:
 * [1] Get the middle-frame between 1 & 2 (m1_2) in the usual way
 *     then the middle-frame between m1_2 & 3 (m1_2_3), etc.
 *     The final result DEPENDS ON the order, i.e., m1_3_2 is different
 *     from m1_2_3, both in position and in orientation!
 * [2] Get the mean z-axis of ALL BASES as the middle-frame z-axis, and
 *     similarly for the x-axis but with a correction for orthogonality.
 *     y-axis follows a right-handed rule. Applying y-first will give
 *     slightly different result; for pan-anti tetrad which is symmetrical,
 *     averaging will give near zero values: ===> big error!
 * So [1] is used all the time ...... */
void pair2mst(long inum_base, long *ivec, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, char **Miscs, double **xyz,
              double **orien, double **org, long **seidx, double *mst_orien,
              double *mst_org, long **htm_water, miscPars *misc_pars,
              FILE *fp) {
  double morg[4], **mst, **xyz_residue;
  long i, ik, j, m, tnum_res, inum = 0;
  long ivec2[BUF512];

  tnum_res = attached_residues(inum_base, ivec, ivec2, seidx, xyz, htm_water,
                               misc_pars);

  mst = dmatrix(1, 3, 1, 3);
  init_dvector(morg, 1, 3, 0.0);
  cehs_average(inum_base, ivec, orien, org, mst, morg);

  if (mst_orien != NULL) { /* return middle frame */
    cpxyz(morg, mst_org);
    mst2orien(mst_orien, 0, mst);
  }

  xyz_residue = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
  for (i = 1; i <= tnum_res; i++) { /* for each residue + attachments */
    ik = ivec2[i];
    for (j = seidx[ik][1]; j <= seidx[ik][2]; j++) {
      m = j - seidx[ik][1] + 1;
      cpxyz(xyz[j], xyz_residue[m]);
    }
    if (!Gvars.ORIGINAL_COORDINATE)
      change_xyz(0, morg, mst, seidx[ik][2] - seidx[ik][1] + 1, xyz_residue);
    pdb_record(seidx[ik][1], seidx[ik][2], &inum, 1, AtomName, ResName, ChainID,
               ResSeq, xyz_residue, Miscs, fp);
  }

  free_dmatrix(mst, 1, 3, 1, 3);
  free_dmatrix(xyz_residue, 1, NUM_RESIDUE_ATOMS, 1, 3);
}

FILE *open_file(char *filename, char *filemode) {
  FILE *fp;

  errno = 0;
  if (filename == NULL)
    filename = "\0";
  if (!strcmp(filename, "stdin"))
    fp = stdin;
  else if (!strcmp(filename, "stdout"))
    fp = stdout;
  else if (!strcmp(filename, "stderr"))
    fp = stderr;
  else {
    fp = fopen(filename, filemode);
    if (fp == NULL) {
      cout << "\n \n \n Filename is: " << filename << " \n \n \n";
      cout << "\n \n \n Filename is: " << filename << " \n \n \n";
      fatal("open_file <%s> failed: %s\n", filename, strerror(errno));
    }
  }
  return fp;
}

long close_file(FILE *fp) {
  long i;

  if (fp == NULL || fp == stdin || fp == stdout || fp == stderr)
    return 0;
  errno = 0;
  i = fclose(fp);
  if (i == EOF)
    fatal("close_file failed: %s\n", strerror(errno));
  return i;
}

long exist_file(char *filename) {
  long status;
  FILE *fp;

  fp = fopen(filename, "r");
  status = (fp != NULL) ? 1 : 0;
  close_file(fp);
  return status;
}

void remove_file(char *filename) {
  if (!exist_file(filename))
    return;

  if (remove(filename))
    fatal("can not remove file: %s\n", filename);
}

void rename_file(char *src, char *dst) {
  if (!exist_file(src))
    fatal("file to be renamed <%s> does NOT exist\n", src);

  if (rename(src, dst))
    fatal("can not rename file <%s> to <%s>\n", src, dst);
}

void copy_file_pointer(FILE *fpi, FILE *fpo, char *msg) {
  char str[BUF512];
  size_t num_bytes;

  while (!feof(fpi)) {
    num_bytes = fread(str, 1, BUF512, fpi);
    if (num_bytes == 0)
      fprintf(stderr, "zero bytes? [%s]\n", msg);
    if (fwrite(str, 1, num_bytes, fpo) != num_bytes)
      fatal("file %s error\n", msg);
  }
}

/* change to upper case, and return string length */
long upperstr(char *a) {
  long nlen = 0;

  while (*a) {
    nlen++;
    if (islower((int)*a))
      *a = toupper((int)*a);
    a++;
  }
  return nlen;
}

/* change to lower case, and return string length */
long lowerstr(char *a) {
  long nlen = 0;

  while (*a) {
    nlen++;
    if (isupper((int)*a))
      *a = tolower((int)*a);
    a++;
  }
  return nlen;
}

char *my_strdup(const char *src) {
  char *dst;

  dst = cvector(0, strlen(src));
  strcpy(dst, src);

  return dst;
}

/* print char 'x' n-times to stream fp */
void print_sep(FILE *fp, char x, long n) {
  long i;

  for (i = 1; i <= n; i++)
    if (fputc(x, fp) == EOF)
      fatal("error writing characters to the stream\n");
  if (fputc('\n', fp) == EOF)
    fatal("error writing '\n' to the stream\n");
}

/* check if '/' exist at the end of the original BDIR */
void check_slash(char *BDIR) {
  char *pchar;
  long n;

  pchar = strrchr(BDIR, '/');
  n = strlen(BDIR);
  if (pchar && (pchar - BDIR != n - 1)) {
    BDIR[n] = '/';
    BDIR[n + 1] = '\0';
  }
}

/* return a pointer 1-char following the last slash or to str w/o '/' */
char *original_name(char *str) {
  char *p_lslash;

  p_lslash = strrchr(str, '/');
  if (p_lslash == NULL)
    return str;
  else
    return p_lslash + 1;
}

long lround(double d) { return (long)((d > 0.0) ? d + 0.5 : d - 0.5); }

/* get rid of the extension in a file name */
void del_extension(char *fullname, char *okname) {
  char *pchar, *bname;
  size_t i;

  bname = original_name(fullname);
  // bname = fullname.substr(fullname.find("/") + 1);

  pchar = strrchr(bname, '.');
  if (pchar == NULL)
    strcpy(okname, bname);
  else {
    i = pchar - bname;
    strncpy(okname, bname, i);
    okname[i] = '\0';
  }
}

/* get the base name without last extension */
void bname_noext(char *src, char *dst) {
  char str[BUF512];

  strcpy(str, original_name(src));
  del_extension(str, dst);
}

/* exit program with an error message */
void fatal(char *fmt, ...) {
  va_list args;

  if (strlen(fmt) > 0) {
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
  }

  exit(1);
}

/* print HEADER, TITLE etc records to parameter output file for easy reference
 */
void print_pdb_title(char *pdbfile, char *chain_list, FILE *fp) {
  char str[BUF512];
  static char *titles[] = {"HEADER", "TITLE ", "AUTHOR", "COMPND",
                           "SOURCE", "EXPDTA", "KEYWDS", "REVDAT",
                           "JRNL  ", "HELIX ", "SHEET ", "TURN  "};
  long i, num, nlen;
  FILE *fpp;

  if (!Gvars.HEADER)
    return;

  num = sizeof titles / sizeof titles[0];

  fpp = open_file(pdbfile, "r");
  while (fgets(str, sizeof str, fpp) != NULL) {
    nlen = upperstr(str);
    if (!strncmp(str, "ATOM  ", 6) || !strncmp(str, "HETATM", 6) ||
        !strncmp(str, "END", 3))
      break;
    if (Gvars.HEADER > true) {
      fprintf(fp, "%s", str);
      continue;
    }

    /* default case */
    if (nlen >= 6) { /* at least 6 characters */
      for (i = 0; i < num - 3; i++)
        if (!strncmp(str, titles[i], 6))
          fprintf(fp, "%s", str);
      if ((!strncmp(str, "HELIX ", 6) || !strncmp(str, "TURN  ", 6)) &&
          (strchr(chain_list, '*') || strchr(chain_list, str[19])))
        fprintf(fp, "%s", str);
      if (!strncmp(str, "SHEET ", 6) &&
          (strchr(chain_list, '*') || strchr(chain_list, str[21])))
        fprintf(fp, "%s", str);
    }
  }
  close_file(fpp);
}

static long is_end_of_structure_to_process(char *str) {
  if (str_pmatch(str, "END")) {
    if (Gvars.ALL_MODEL)
      return str_pmatch(str, "ENDMDL") ? false : true;
    else /* also matches ENDMDL */
      return true;
  } else
    return false;
}

static double get_occupancy(long nlen, char *str, char *pdbfile) {
  char temp[BUF512];
  double occupancy;

  if (nlen < 60 || !Gvars.OCCUPANCY) /* no checking for occupancy */
    return 1.0;

  strncpy(temp, str + 54, 6); /* occupancy */
  temp[6] = '\0';
  if (sscanf(temp, "%6lf", &occupancy) != 1)
    fatal("error reading occupancy in file [%s]: '%s'\n", pdbfile, str);

  return occupancy;
}

/* The last column of z-coordinate is 54. However, some non-compliant
 * PDB could be less than 54 column per ATOM/HETATM record. So here
 * use 52 to account for such special cases */
#define Zcol 52

/* number of atom records in a PDB file */
long number_of_atoms(char *pdbfile, long hetatm, char *ALT_LIST) {
  char str[BUF512], *pchar;
  long n = 0, nlen;
  FILE *fp;

  fp = open_file(pdbfile, "r");
  while (fgets(str, sizeof str, fp) != NULL) {
    if ((pchar = strrchr(str, '\n')) != NULL)
      str[pchar - str] = '\0';
    nlen = upperstr(str);

    if (is_end_of_structure_to_process(str))
      break;

    if (nlen >= Zcol &&
        (!strncmp(str, "ATOM", 4) || (hetatm && !strncmp(str, "HETATM", 6))) &&
        (strchr(ALT_LIST, '*') || strchr(ALT_LIST, str[16])) &&
        get_occupancy(nlen, str, pdbfile) > 0)
      n++;
  }
  close_file(fp);

  return n;
}

/* read in a PDB file and do some processing
 * Miscs[][NMISC]: H/A, altLoc, iCode, occ./tempFac./segID/element/charge
 *           col#   0     1       2       3-28 [combined together]
 * 29-32, 4 characters for Model Numbers (in PDB MODEL record, col: 11-14)
 * Element symbol: column# 25-26, right justified
 * 80 - 54 = 26 characters from ATOM/HETATM records 26 + 2 = 28
 *      set the 29th column to '\0' for end of the ATOM/HETATM record
 * 30 -- 33 for the 4-character model number; 34 is set to '\0'
 *   so we have NMISC = 34 */
long read_pdb(char *pdbfile, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs,
              long hetatm, char *ALT_LIST) {
  char *p0, *pn;
  char rname[4], str[BUF512], str0[BUF512], temp[BUF512];
  double occupancy;
  long i, n = 0, nlen, k, occ_chk = false, modelNum = 0;
  FILE *fp;

  fp = open_file(pdbfile, "r");
  while ((p0 = my_getline(fp)) != NULL) {
    strcpy(str, p0);
    free(p0);

    strcpy(str0, str);
    nlen = upperstr(str);
    if (!strncmp(str, "MODEL ", 6) /* model record */
        && (sscanf(str, "%*s %ld", &k) == 1))
      modelNum = k;

    if (is_end_of_structure_to_process(str))
      break;

    if (nlen >= Zcol && /* same as in number_of_atoms */
        (!strncmp(str, "ATOM", 4) || (hetatm && !strncmp(str, "HETATM", 6))) &&
        (strchr(ALT_LIST, '*') || strchr(ALT_LIST, str[16]))) {

      occupancy = get_occupancy(nlen, str, pdbfile);
      if (occupancy <= 0) { /* ignore 0-occupancy atom */
        if (!occ_chk) {
          occ_chk = true;
          fprintf(stderr, "[i] File '%s' with atom occupancy <= 0 [%s]\n",
                  pdbfile, str);
        }
        continue;
      }

      n++;

      if (AtomSNum != NULL) { /* read in original atom serial number */
        strncpy(temp, str + 6, 5);
        temp[5] = '\0'; /* 5 digits */
        if (sscanf(temp, "%5ld", &AtomSNum[n]) != 1) {
          fprintf(stderr, "Atom serial number ? %s ?\n", str);
          AtomSNum[n] = n; /* sequential number */
        }
      }

      strncpy(AtomName[n], str + 12, 4);
      AtomName[n][4] = '\0';

      if (Gvars.AtomName0 && Gvars.Name0)
        strcpy(Gvars.AtomName0[n], AtomName[n]);

      /* for very special case 'P   ', 'N9  ' etc. Fei XU: also "MG  " 355d */
      pn = AtomName[n]; /* as a shorthand */
      if ((pn[0] != ' ' && !isdigit((int)pn[0])) &&
          (pn[1] == ' ' || isdigit((int)pn[1])) && pn[3] == ' ') {
        strncpy(temp, pn, 3);
        temp[3] = '\0';
        sprintf(pn, " %s", temp);
      } else if (pn[0] == ' ' && pn[1] == ' ' &&
                 isdigit((int)pn[3])) { /* as in "  N1" */
        strncpy(temp, pn + 2, 2);
        temp[2] = '\0';
        sprintf(pn, " %s ", temp);
      } else if (is_equal_string(pn, "   P") || is_equal_string(pn, "P   "))
        strcpy(pn, " P  ");
      else if (is_equal_string(pn, "OP1 "))
        strcpy(pn, " O1P");
      else if (is_equal_string(pn, "OP2 "))
        strcpy(pn, " O2P");

      if (AtomName[n][3] == '*') /* * to ' */
        AtomName[n][3] = '\'';
      if (!strcmp(AtomName[n], " O1'")) /* O1' to O4' */
        strcpy(AtomName[n], " O4'");
      if (!strcmp(AtomName[n], " OL ") ||
          !strcmp(AtomName[n], " OP1")) /* OL/OP1 to O1P */
        strcpy(AtomName[n], " O1P");
      if (!strcmp(AtomName[n], " OR ") ||
          !strcmp(AtomName[n], " OP2")) /* OR/OP2 to O2P */
        strcpy(AtomName[n], " O2P");
      if (!strcmp(AtomName[n], " OP3")) /* OP3 to O3P */
        strcpy(AtomName[n], " O3P");
      if (!strcmp(AtomName[n], " C5A")) /* C5A to C5M */
        strcpy(AtomName[n], " C5M");
      if (!strcmp(AtomName[n], " O5T")) /* terminal O5' */
        strcpy(AtomName[n], " O5'");
      if (!strcmp(AtomName[n], " O3T")) /* terminal O3' */
        strcpy(AtomName[n], " O3'");

      strncpy(rname, str + 17, 3); /* residue name */
      rname[3] = '\0';
      if (Gvars.ResName0 && Gvars.Name0)
        strcpy(Gvars.ResName0[n], rname);

      /* delete ending spaces as in "C  " */
      for (i = 1; i <= 2; i++)
        if (rname[2] == ' ') {
          rname[2] = rname[1];
          rname[1] = rname[0];
          rname[0] = ' ';
        }
      if (rname[2] == ' ')
        fatal("==> residue name [%s] field empty <==\n", str);
      strcpy(ResName[n], rname);

      ChainID[n] = Gvars.CHAIN_CASE ? str0[21] : str[21];
      strncpy(temp, str + 22, 4);
      temp[4] = '\0';
      if (sscanf(temp, "%4ld", &ResSeq[n]) != 1) {
        fprintf(stderr, "residue #? ==> %.54s\n", str);
        ResSeq[n] = 9999;
      }
      strncpy(temp, str + 30, 8);
      temp[8] = '\0';
      if (sscanf(temp, "%8lf", &xyz[n][1]) != 1)
        fatal("error reading x-coordinate\n");
      strncpy(temp, str + 38, 8);
      temp[8] = '\0';
      if (sscanf(temp, "%8lf", &xyz[n][2]) != 1)
        fatal("error reading y-coordinate\n");
      strncpy(temp, str + 46, 8);
      temp[8] = '\0';
      if (sscanf(temp, "%8lf", &xyz[n][3]) != 1)
        fatal("error reading z-coordinate\n");
      if (Miscs != NULL) {
        Miscs[n][0] = str[0];  /* H for HETATM, A for ATOM */
        Miscs[n][1] = str[16]; /* alternative location indicator */
        Miscs[n][2] = str[26]; /* code of insertion residues */
        if (nlen >= 54)
          strncpy(Miscs[n] + 3, str + 54, 26); /* upto the 80th column */
        Miscs[n][29] = '\0';
        sprintf(Miscs[n] + 30, "%4.4ld", modelNum);
        Miscs[n][NMISC] = '\0'; /* just to make sure */
      }
    }
  }
  close_file(fp);

  if (!n)
    fprintf(stderr, "PDB file <%s> has NO ATOM/HETATM records\n", pdbfile);

  return n;
}

/* free all PDB relevant arrays to make code clean and short */
void free_pdb(long num, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs) {
  if (AtomSNum != NULL)
    free_lvector(AtomSNum, 1, num);
  if (AtomName != NULL)
    free_cmatrix(AtomName, 1, num, 0, 4);
  if (ResName != NULL)
    free_cmatrix(ResName, 1, num, 0, 3);
  if (ChainID != NULL)
    free_cvector(ChainID, 1, num);
  if (ResSeq != NULL)
    free_lvector(ResSeq, 1, num);
  if (xyz != NULL)
    free_dmatrix(xyz, 1, num, 1, 3);
  if (Miscs != NULL)
    free_cmatrix(Miscs, 1, num, 0, NMISC);
}

/* reset xyz coordinates for PDB and ALCHEMY formats */
void reset_xyz(long num, double **xyz, char *fmt) {
  double ave_xyz[4], max_xyz[4], min_xyz[4];

  /* check the range of the coordinates */
  max_dmatrix(xyz, num, 3, max_xyz);
  min_dmatrix(xyz, num, 3, min_xyz);
  ave_dmatrix(xyz, num, 3, ave_xyz);

  if (max_dvector(max_xyz, 1, 3) > 9999.99) {
    fprintf(stderr,
            "xyz coordinate over %s limit. reset origin"
            " to geometrical center\n",
            fmt);
    move_position(xyz, num, 3, ave_xyz);
  } else if (min_dvector(min_xyz, 1, 3) < -999.99) {
    fprintf(stderr,
            "xyz coordinate under %s limit. reset origin"
            " to minimum xyz coordinates\n",
            fmt);
    move_position(xyz, num, 3, min_xyz);
  }
}

void deduce_misc(char **Miscs, char **AtomName, long i, char *str) {
  static char asym[3];

  if (Miscs == NULL || is_equal_string(Miscs[i], "A  ")) {
    aname2asym(AtomName[i], asym, Gvars.NUM_SATOM, Gvars.ATOMLIST);
    /* ATOM/HETATM, altLoc, iCode, occupancy, temp-factor, unused 67..76,
     * atom-symbol, and charge */
    sprintf(str, "A  %6.2f%6.2f          %2.2s  ", 1.0, 1.0, asym);
  } else
    strcpy(str, Miscs[i]);
}

static void cvt_3letter_nts(char *rname) {
  if (is_equal_string(rname, "  A") || is_equal_string(rname, " DA"))
    strcpy(rname, "ADE");
  else if (is_equal_string(rname, "  C") || is_equal_string(rname, " DC"))
    strcpy(rname, "CYT");
  else if (is_equal_string(rname, "  G") || is_equal_string(rname, " DG"))
    strcpy(rname, "GUA");
  else if (is_equal_string(rname, "  T") || is_equal_string(rname, " DT"))
    strcpy(rname, "THY");
  else if (is_equal_string(rname, "  U"))
    strcpy(rname, "URA");
}

long is_dna_with_backbone(long ib, long ie, char **AtomName) {
  long i, P = false;

  for (i = ib; i <= ie; i++) {
    if (is_equal_string(AtomName[i], " O2'"))
      return false; /* taken as RNA */
    if (!P && is_equal_string(AtomName[i], " P  "))
      P = true;
  }

  return P;
}

static void cvt_pdbv3_name(char *aname) {
  if (is_equal_string(aname, " O1P"))
    strcpy(aname, " OP1");
  else if (is_equal_string(aname, " O2P"))
    strcpy(aname, " OP2");
  else if (is_equal_string(aname, " C5M"))
    strcpy(aname, " C7 ");
}

static void cvt_pdbv3_dna(char *rname) {
  if (is_equal_string(rname, "  A"))
    strcpy(rname, " DA");
  else if (is_equal_string(rname, "  C"))
    strcpy(rname, " DC");
  else if (is_equal_string(rname, "  G"))
    strcpy(rname, " DG");
  else if (is_equal_string(rname, "  T"))
    strcpy(rname, " DT");
}

void normalize_resName_atomName(long is_dna, const char *rname0,
                                const char *aname0, char *rname, char *aname) {
  strcpy(aname, aname0);
  strcpy(rname, rname0);

  if (Gvars.PDBV3) {
    cvt_pdbv3_name(aname);
    if (is_dna)
      cvt_pdbv3_dna(rname);

  } else if (Gvars.THREE_LETTER_NTS)
    cvt_3letter_nts(rname);
}

/* write out ATOM and HETATM record: xyz could be 1 to [ie - ib + 1] */
void pdb_record(long ib, long ie, long *inum, long idx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, double **xyz,
                char **Miscs, FILE *fp) {
  char rname[4], aname[5], str[BUF512];
  long i, j, is_dna;

  is_dna = Gvars.PDBV3 && is_dna_with_backbone(ib, ie, AtomName);

  for (i = ib; i <= ie; i++) {
    deduce_misc(Miscs, AtomName, i, str);
    if (Gvars.AtomName0) {
      strcpy(aname, Gvars.AtomName0[i]);
      strcpy(rname, Gvars.ResName0[i]);
    } else
      normalize_resName_atomName(is_dna, ResName[i], AtomName[i], rname, aname);

    j = (idx) ? i - ib + 1 : i;
    fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
            (str[0] == 'A') ? "ATOM  " : "HETATM", ++*inum, aname, str[1],
            rname, ChainID[i], ResSeq[i], str[2], xyz[j][1], xyz[j][2],
            xyz[j][3], str + 3);
  }
}

static void write_cif_header(FILE *fp) {
  fprintf(fp, "data_x3dna\n");
  fprintf(fp, "# mmCIF output file generated by 3DNA (xiangjun@x3dna.org)\n");

  fprintf(fp, "loop_\n");

  fprintf(fp, "_atom_site.group_PDB\n"     /* ATOM/HETATM */
              "_atom_site.id\n"            /* serial number */
              "_atom_site.type_symbol\n"); /* atom symbol */

  fprintf(fp, "_atom_site.label_atom_id\n"  /* atom name */
              "_atom_site.label_alt_id\n"   /* alternative location */
              "_atom_site.label_comp_id\n"  /* residue name */
              "_atom_site.label_asym_id\n"  /* chain id */
              "_atom_site.label_seq_id\n"); /* sequence number, renumber from 1
                                               for each chain */

  fprintf(fp, "_atom_site.pdbx_PDB_ins_code\n" /* insertion code */
              "_atom_site.Cartn_x\n"           /* x-coordinate */
              "_atom_site.Cartn_y\n"           /* y-coordinate */
              "_atom_site.Cartn_z\n");         /* z-coordinate */

  fprintf(fp, "_atom_site.occupancy\n"            /* occupancy */
              "_atom_site.B_iso_or_equiv\n"       /* B-factor */
              "_atom_site.pdbx_formal_charge\n"); /* formal charge */

  fprintf(fp, "_atom_site.auth_seq_id\n"    /* sequence number */
              "_atom_site.auth_comp_id\n"   /* residue name */
              "_atom_site.auth_asym_id\n"   /* chain id */
              "_atom_site.auth_atom_id\n"); /* atom name */

  fprintf(fp, "_atom_site.pdbx_PDB_model_num\n"); /* model number */
}

void move_position(double **d, long nr, long nc, double *mpos) {
  long i, j;

  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
      d[i][j] -= mpos[j];
}

/* number of residues, and starting-ending indexes for each */
long **residue_idx(long num, long *ResSeq, char **Miscs, char *ChainID,
                   char **ResName, long *num_residue) {
  char iCode;
  char **bidx;
  long i, n, **seidx, *temp;

  bidx = cmatrix(1, num, 0, 12); /* normally 9 */
  temp = lvector(1, num);

  for (i = 1; i <= num; i++) {
    iCode = (Miscs == NULL) ? ' ' : Miscs[i][2];
    sprintf(bidx[i], "%3s%c%4ld%c", ResName[i], ChainID[i], ResSeq[i], iCode);
  }
  for (i = 1; i < num; i++)
    temp[i] = strcmp(bidx[i + 1], bidx[i]) ? 1 : 0;
  temp[num] = 1;

  n = 0; /* get number of residues */
  for (i = 1; i <= num; i++)
    if (temp[i])
      ++n;

  seidx = lmatrix(1, n, 1, 2); /* allocate spaces */
  n = 0;
  for (i = 1; i <= num; i++)
    if (temp[i])
      seidx[++n][2] = i;
  for (i = 2; i <= n; i++)
    seidx[i][1] = seidx[i - 1][2] + 1;
  seidx[1][1] = 1;

  *num_residue = n;

  free_cmatrix(bidx, 1, num, 0, 12);
  free_lvector(temp, 1, num);

  return seidx;
}

static void set_U_C5M(char **AtomName, double **xyz, long ib, long ie) {
  long C5, C5M, C7;

  C5 = find_1st_atom(" C5 ", AtomName, ib, ie, "");
  C5M = find_1st_atom(" C5M", AtomName, ib, ie, "");
  C7 = find_1st_atom(" C7 ", AtomName, ib, ie, "");

  if (C5 && C7 && !C5M && (p1p2_dist(xyz[C5], xyz[C7]) < 2.0))
    strcpy(AtomName[C7], " C5M");
}

void atom_metal(long num_atoms, char **AtomName, long *is_metal) {
  static char *metals[] = {
      "LI", "BE", "NA", "MG", "AL", " K", "CA", "SC", "TI", " V", "CR", "MN",
      "FE", "CO", "NI", "CU", "ZN", "GA", "RB", "SR", " Y", "ZR", "NB", "MO",
      "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "CS", "BA", "HF", "TA",
      " W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "FR", "RA",
      "RF", "DB", "SG", "BH", "HS", "MT", "LA", "CE", "PR", "ND", "PM", "SM",
      "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "AC", "TH", "PA",
      " U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR"};

#if 0
    static char *non_metals[] = {
        " H", "HE",
        " B", " C", " N", " O", " F", "NE",
        "SI", " P", " S", "CL", "AR",
        "GE", "AS", "SE", "BR", "KR",
        "SB", "TE", " I", "XE",
        "PO", "AT", "RN"
    };
#endif

  char atom_sym[3];
  long i, num_metal;

  num_metal = sizeof metals / sizeof metals[0] - 1; /* minus - 1 */

  for (i = 1; i <= num_atoms; i++) {
    aname2asym(AtomName[i], atom_sym, Gvars.NUM_SATOM, Gvars.ATOMLIST);
    is_metal[i] = (num_strmatch(atom_sym, metals, 0, num_metal)) ? true : false;
  }
}

/* 2o8b_C30F has RMSD 0.24; normally < 0.1 */
static double check_nt_type_by_rmsd(long *idx, long C1_prime, double **xyz) {
  /* idx[] RA_LIST: " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8
   * ", " N9 " */
  static double xyz_ring[][4] = {
      /* based on Atomic_G.pdb */
      {XBIG, -1.265, 3.177, 0.000}, /* C4 -- 0 */
      {XBIG, -2.342, 2.364, 0.001}, /* N3 -- 1 */
      {XBIG, -1.999, 1.087, 0.000}, /* C2 -- 2 */
      {XBIG, -0.700, 0.641, 0.000}, /* N1 -- 3 */
      {XBIG, 0.424, 1.460, 0.000},  /* C6 -- 4 */
      {XBIG, 0.071, 2.833, 0.000},  /* C5 -- 5 */
      {XBIG, 0.870, 3.969, 0.000},  /* N7 -- 6 */
      {XBIG, 0.023, 4.962, 0.000},  /* C8 -- 7 */
      {XBIG, -1.289, 4.551, 0.000}  /* N9 -- 8 */
  };
  double **xyz1, **xyz2, **fitted_xyz, **R;
  double rmsd, org[4];
  long i, k = 0, nN = 0, num = 9;

  /* large enough to hold everything: 3 * 9 + 3 = 30 */
  xyz1 = dmatrix(1, BUF32, 1, 3);
  xyz2 = &xyz1[num + 1];
  fitted_xyz = &xyz1[2 * num + 1];
  R = &xyz1[3 * num + 1];

  for (i = 0; i < num; i++) {
    if (!idx[i])
      continue;
    k++;
    if (i == 1 || i == 3 || i == 6 || i == 8)
      nN++;

    cpxyz(xyz[idx[i]], xyz1[k]);
    cpxyz(xyz_ring[i], xyz2[k]);
  }

  if (!nN && !C1_prime)
    rmsd = DUMMY;
  else
    rmsd = ls_fitting(xyz1, xyz2, k, fitted_xyz, R, org);

  free_dmatrix(xyz1, 1, DUMMY, 1, DUMMY);

  return rmsd;
}

/* identifying a residue as follows:
 *  R-base  Y-base  amino-acid, others [default]
 *   +1        0        -1        -2 [default] */
long residue_ident(char **AtomName, double **xyz, char **Miscs, long ib,
                   long ie) {
  static char *RingAtom[] = {RA_LIST};
  double rmsd = XBIG, dcrt = 2.0;
  long i, n, k = 0, kr = 0, num = 9, idx[16] = {0};
  long CA, C, C1_prime;

  for (i = 0; i < num; i++) { /* nine atoms */
    n = find_1st_atom(RingAtom[i], AtomName, ib, ie, "");
    if (n) {
      k++;
      if (i >= 6) /* N7/C8/N9 */
        kr++;
      idx[i] = n;
    } else
      idx[i] = 0; /* make explicit assignment */
  }

  C1_prime = find_1st_atom(" C1'", AtomName, ib, ie, "");
  if (k >= 3) /* at least 3 ring atoms */
    rmsd = check_nt_type_by_rmsd(idx, C1_prime, xyz);

  if (rmsd != DUMMY && rmsd <= Gvars.NT_CUTOFF) {
    if (kr)
      return 1;
    set_U_C5M(AtomName, xyz, ib, ie);
    return 0;

  } else if (kr) { /* added on 2019-02-18 */
    kr = 0;
    for (i = 6; i < num; i++)
      idx[i] = 0;
    k = 0;
    for (i = 0; i < 6; i++)
      if (idx[i])
        k++;
    if (k >= 3)
      rmsd = check_nt_type_by_rmsd(idx, C1_prime, xyz);
    if (rmsd != DUMMY && rmsd <= Gvars.NT_CUTOFF) {
      set_U_C5M(AtomName, xyz, ib, ie);
      return 0;
    }
  }

  CA = find_1st_atom(" CA ", AtomName, ib, ie, "");
  C = find_1st_atom(" C  ", AtomName, ib, ie, "");
  if (!C) /* if C does not exist, use N */
    C = find_1st_atom(" N  ", AtomName, ib, ie, "");
  if (CA && C && Miscs[CA][0] == 'A' && Miscs[C][0] == 'A' && /* ATOM record */
      within_limits(xyz[CA], xyz[C], 0, dcrt))
    return -1;

  return -2;
}

void normalize_atom_symbol(char *asym) {
  upperstr(asym);

  if (strlen(asym) == 1) { /* "H" --> " H" */
    asym[1] = asym[0];
    asym[0] = ' ';
    asym[2] = '\0';
  }

  if (is_equal_string(asym, " D"))
    strcpy(asym, " H");
}

/* get the correspondence between 4-letter atom name to atomic symbol */
void get_atomlist(char **atomlist, long *num_sa) {
  char BDIR[BUF1K], str[BUF512], aname4[BUF512], asym[BUF512];
  long n = 0;
  FILE *fp;

  get_BDIR(BDIR, ATOM_FILE);
  strcat(BDIR, ATOM_FILE);
  if (Gvars.VERBOSE)
    fprintf(stderr, " ...... reading file: %s ...... \n", ATOM_FILE);

  fp = open_file(BDIR, "r");

  while ((fgets(str, sizeof str, fp) != NULL)) {
    if (str[0] == '#')
      continue;
    if (sscanf(str, "%s %s", aname4, asym) != 2)
      continue;
    if (aname4[0] == '#' || asym[0] == '#')
      continue;

    if (strlen(aname4) != 4) {
      if (Gvars.VERBOSE)
        fprintf(stderr,
                "atom name must be 4-char long with only [.A-Z]: <%s>\n",
                aname4);
      continue;
    }

    if (strlen(asym) != 1 && strlen(asym) != 2) {
      if (Gvars.VERBOSE)
        fprintf(stderr,
                "atomic symbol must be 1/2-char with only [A-Z]: <%s>\n", asym);
      continue;
    }

    normalize_atom_symbol(asym);
    if (!num_strmatch(asym, Gvars.ATOM_NAMES, 0, Gvars.NUM_ELE)) {
      if (Gvars.VERBOSE)
        fprintf(stderr, "skip invalid atom symbol <%s : %s>\n", asym, aname4);
      continue;
    }

    upperstr(aname4);

    if (++n > BUFBIG)
      fatal("too many atom types\n");

    sprintf(atomlist[n], "%4.4s%2.2s", aname4, asym);
  }

  close_file(fp);

  *num_sa = n;
}

long has_atom_name(long ib, long ie, char **AtomName, char *aname) {
  long i;

  for (i = ib; i <= ie; i++)
    if (strcmp(AtomName[i], aname) == 0)
      return 1L;

  return 0L;
}

/* get one letter base name according to 3-letter to 1-letter conversion */
static char base_ident(long ib, long ie, char **AtomName, double **xyz,
                       long isR, char *rname, char *idmsg, long num_sb,
                       char **baselist) {
  char c = 'X';
  long i;

  for (i = 1; i <= num_sb; i++)
    if (!strncmp(rname, baselist[i], 3))
      break;

  if (i > num_sb) {
    if (isR) {
      if (has_atom_name(ib, ie, AtomName,
                        " O6 ")) /* instead of N2: 2013oct13/14 */
        c = 'g';
      else if (!has_atom_name(ib, ie, AtomName, " N6 ") &&
               has_atom_name(ib, ie, AtomName, " N2 "))
        c = 'g';
      else
        c = 'a';
    } else {
      if (has_atom_name(ib, ie, AtomName, " N4 "))
        c = 'c';
      else if (has_atom_name(ib, ie, AtomName, " C5M"))
        c = 't';
      else
        c = 'u';
      {
        long c1p, n1, c5;
        c1p = find_1st_atom(" C1'", AtomName, ib, ie, "");
        n1 = find_1st_atom(" N1 ", AtomName, ib, ie, "");
        c5 = find_1st_atom(" C5 ", AtomName, ib, ie, "");
        if (c1p && n1 && c5 && p1p2_dist(xyz[c1p], xyz[n1]) > 2.0 &&
            p1p2_dist(xyz[c1p], xyz[c5]) <= 2.0)
          c = 'p'; /* as in 3TD of 5afi */
      }
    }
    fprintf(stderr, "Match '%s' to '%c' for %s\n", rname, c, idmsg);
    fprintf(stderr,
            "    check it & consider to add line '%s     %c' to file <%s>\n",
            rname, c, BASE_FILE);
  } else {
    c = baselist[i][3];
    if (!isupper((int)c) || (c == 'P'))
      fprintf(stderr, "[i] uncommon %s assigned to: %c\n", idmsg, c);
  }

  return c;
}

/* get the correspondence between 3-letter and 1-letter base names */
void get_baselist(char **baselist, long *num_sb) {
  char BDIR[BUF1K], str[BUF512], base3[BUF512], base1[BUF512];
  long n = 0;
  FILE *fp;

  get_BDIR(BDIR, BASE_FILE);
  strcat(BDIR, BASE_FILE);
  if (Gvars.VERBOSE)
    fprintf(stderr, " ...... reading file: %s ...... \n", BASE_FILE);

  fp = open_file(BDIR, "r");

  while ((fgets(str, sizeof str, fp) != NULL)) {
    if (str[0] == '#')
      continue;
    if (sscanf(str, "%s %s", base3, base1) != 2)
      continue;
    if (base3[0] == '#' || base1[0] == '#')
      continue;

    if (strlen(base3) > 3 || strlen(base1) != 1) {
      if (Gvars.VERBOSE)
        fprintf(stderr, "ignoring unacceptable format: %s\n", str);
      continue;
    }

    if (++n > BUFBIG)
      fatal("too many base types\n");

    upperstr(base3);
    sprintf(baselist[n], "%3.3s%c", base3, base1[0]);
  }

  close_file(fp);

  *num_sb = n;
}

/* get base sequence, similar to GET_BPSEQ */
void get_seq(long num_residue, long **seidx, char **AtomName, char **ResName,
             char *ChainID, long *ResSeq, char **Miscs, double **xyz,
             char *bseq, long *RY) {
  char idmsg[BUF512];
  long i, ib, ie;

  for (i = 1; i <= num_residue; i++) {
    ib = seidx[i][1];
    ie = seidx[i][2];
    RY[i] = residue_ident(AtomName, xyz, Miscs, ib, ie);
    if (RY[i] >= 0) {
      sprintf(idmsg, "residue %3s %4ld%c on chain %c [#%ld]", ResName[ib],
              ResSeq[ib], Miscs[ib][2], ChainID[ib], i);
      bseq[i] = base_ident(ib, ie, AtomName, xyz, RY[i], ResName[ib], idmsg,
                           Gvars.NUM_SBASE, Gvars.BASELIST);
    }
  }
}

/* get base (pair) sequence, change (+) modified residue to lower case */
void get_bpseq(long ds, long num_bp, long **pair_num, long **seidx,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, char **bp_seq, long *RY) {
  char idmsg[BUF512];
  long i, ib, ie, j, rnum;

  for (i = 1; i <= ds; i++) {
    for (j = 1; j <= num_bp; j++) {
      rnum = pair_num[i][j];
      ib = seidx[rnum][1];
      ie = seidx[rnum][2];
      RY[rnum] = residue_ident(AtomName, xyz, Miscs, ib, ie);
      sprintf(idmsg, "residue %3s %4ld%c on chain %c [#%ld]", ResName[ib],
              ResSeq[ib], Miscs[ib][2], ChainID[ib], rnum);
      if (RY[rnum] >= 0)
        bp_seq[i][j] = base_ident(ib, ie, AtomName, xyz, RY[rnum], ResName[ib],
                                  idmsg, Gvars.NUM_SBASE, Gvars.BASELIST);
      else
        fatal("Non-base: %s\n", idmsg);
    }
  }
}

/*  return number of matchs of str in strmat */
long num_strmatch(char *str, char **strmat, long nb, long ne) {
  long i, num = 0;

  for (i = nb; i <= ne; i++)
    if (!strcmp(str, strmat[i]))
      num++;

  return num;
}

void get_idmsg(char *rname, char cid, long snum, char icode, char *idmsg) {
  sprintf(idmsg, ": residue name '%s', chain %c, number [%4ld%c]", rname, cid,
          snum, icode);
}

/* return index of the first match, or 0 for no-match */
long find_1st_atom(char *str, char **strmat, long nb, long ne, char *idmsg) {
  char aname[BUF32];
  long i, num;

  num = num_strmatch(str, strmat, nb, ne);

  if (!num) {
    if (strcmp(idmsg, "")) {
      strcpy(aname, str);
      if (Gvars.PDBV3) {
        if (is_equal_string(" O1P", str))
          strcpy(aname, " OP1");
        else if (is_equal_string(" O2P", str))
          strcpy(aname, " OP2");
        else if (is_equal_string(" C5M", str))
          strcpy(aname, " C7 ");
      }
      fprintf(stderr, "[i] missing '%s' atom %s\n", aname, idmsg);
    }
    return 0;
  }

  if (num > 1 && strcmp(idmsg, "")) {
    fprintf(stderr, "more than one %s atoms %s\n", str, idmsg);
    fprintf(stderr, "   *****the first atom is used*****\n");
  }

  for (i = nb; i <= ne; i++)
    if (!strcmp(str, strmat[i]))
      break;

  return i;
}

/* get torsion angle a-b-c-d in degrees */
double torsion(double **d) {
  double ang_deg, **vec3;
  long i;

  vec3 = dmatrix(1, 3, 1, 3);

  for (i = 1; i <= 3; i++) {
    ddxyz(d[i], d[i + 1], vec3[i]);
    if (veclen(vec3[i]) > BOND_UPPER_LIMIT) {
      ang_deg = EMPTY_NUMBER; /* not directly linked */
      goto RTN_HERE;
    }
  }
  negate_xyz(vec3[1]); /* b-->a */
  ang_deg = vec_ang(vec3[1], vec3[3], vec3[2]);

RTN_HERE:
  free_dmatrix(vec3, 1, 3, 1, 3);

  return ang_deg;
}

/* get torsion angle a-b-c-d in degrees: no consideration of breaks */
double torsion2(double **d) {
  double ang_deg, **vec3;
  long i;

  vec3 = dmatrix(1, 3, 1, 3);
  for (i = 1; i <= 3; i++)
    ddxyz(d[i], d[i + 1], vec3[i]);
  negate_xyz(vec3[1]); /* b-->a */
  ang_deg = vec_ang(vec3[1], vec3[3], vec3[2]);
  free_dmatrix(vec3, 1, 3, 1, 3);

  return ang_deg;
}

/* search the directory containing standard bases & parameter files:
   (1) current directory
   (2) directory defined by the environment variable "X3DNA" */
void get_BDIR(char *BDIR, char *filename) {
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp != NULL)
    strcpy(BDIR, "./"); /* current directory */
  else
    sprintf(BDIR, "%sconfig/", Gvars.X3DNA_HOMEDIR);
  close_file(fp);

  if (Gvars.VERBOSE)
    fprintf(stderr, "\n ...... %s%s...... \n", BDIR, filename);
}

/* transform xyz coordinates so that "haxis" is aligned to the z-axis */
void align2zaxis(long num, double *haxis, double **rotmat, double **xyz,
                 double **xyzH) {
  double z[4] = {EMPTY_NUMBER, 0.0, 0.0, 1.0};
  double ang_deg, hinge[4], **rotmatT;

  rotmatT = dmatrix(1, 3, 1, 3);

  cross(haxis, z, hinge);
  ang_deg = magang(haxis, z);
  arb_rotation(hinge, ang_deg, rotmat);
  transpose_matrix(rotmat, 3, 3, rotmatT);
  multi_matrix(xyz, num, 3, rotmatT, 3, 3, xyzH);

  free_dmatrix(rotmatT, 1, 3, 1, 3);
}

/* calculate the covariance matrix between two matrices */
void cov_matrix(double **a, double **b, long nr, long nc, double **cmtx) {
  double ave_a[4], ave_b[4];
  double **ta, **ta_x_b;
  long i, j;

  ave_dmatrix(a, nr, nc, ave_a);
  ave_dmatrix(b, nr, nc, ave_b);

  ta = dmatrix(1, nc, 1, nr);     /* transpose of a */
  ta_x_b = dmatrix(1, nc, 1, nc); /* transpose-a multiply b */

  transpose_matrix(a, nr, nc, ta);
  multi_matrix(ta, nc, nr, b, nr, nc, ta_x_b);

  for (i = 1; i <= nc; i++)
    for (j = 1; j <= nc; j++)
      cmtx[i][j] = (ta_x_b[i][j] - ave_a[i] * ave_b[j] * nr) / (nr - 1);

  free_dmatrix(ta, 1, nc, 1, nr);
  free_dmatrix(ta_x_b, 1, nc, 1, nc);
}

/* least-squares fitting between two structures */
double ls_fitting(double **sxyz, double **exyz, long n, double **fitted_xyz,
                  double **R, double *orgi) {
  double temp, rms_value;
  double ave_exyz[4], ave_sxyz[4], D[5];
  double **N, **U, **V;
  long i, j;

  if (n < 3)
    fatal("too few atoms for least-squares fitting\n");

  /* get the covariance matrix U */
  U = dmatrix(1, 3, 1, 3);
  cov_matrix(sxyz, exyz, n, 3, U);

  /* get 4-by-4 symmetric matrix N */
  N = dmatrix(1, 4, 1, 4);
  N[1][1] = U[1][1] + U[2][2] + U[3][3];
  N[2][2] = U[1][1] - U[2][2] - U[3][3];
  N[3][3] = -U[1][1] + U[2][2] - U[3][3];
  N[4][4] = -U[1][1] - U[2][2] + U[3][3];
  N[1][2] = U[2][3] - U[3][2];
  N[2][1] = N[1][2];
  N[1][3] = U[3][1] - U[1][3];
  N[3][1] = N[1][3];
  N[1][4] = U[1][2] - U[2][1];
  N[4][1] = N[1][4];
  N[2][3] = U[1][2] + U[2][1];
  N[3][2] = N[2][3];
  N[2][4] = U[3][1] + U[1][3];
  N[4][2] = N[2][4];
  N[3][4] = U[2][3] + U[3][2];
  N[4][3] = N[3][4];

  /* get N's eigenvalues and eigenvectors */
  V = dmatrix(1, 4, 1, 4);
  jacobi(N, 4, D, V);

  /* get the rotation matrix */
  for (i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      N[i][j] = V[i][4] * V[j][4];
  R[1][1] = N[1][1] + N[2][2] - N[3][3] - N[4][4];
  R[1][2] = 2 * (N[2][3] - N[1][4]);
  R[1][3] = 2 * (N[2][4] + N[1][3]);
  R[2][1] = 2 * (N[3][2] + N[1][4]);
  R[2][2] = N[1][1] - N[2][2] + N[3][3] - N[4][4];
  R[2][3] = 2 * (N[3][4] - N[1][2]);
  R[3][1] = 2 * (N[4][2] - N[1][3]);
  R[3][2] = 2 * (N[4][3] + N[1][2]);
  R[3][3] = N[1][1] - N[2][2] - N[3][3] + N[4][4];

  ave_dmatrix(sxyz, n, 3, ave_sxyz);
  ave_dmatrix(exyz, n, 3, ave_exyz);

  /* fitted sxyz origin */
  for (i = 1; i <= 3; i++)
    orgi[i] = ave_exyz[i] - dot(ave_sxyz, R[i]);

  /* fitted sxyz coordinates */
  for (i = 1; i <= n; i++)
    for (j = 1; j <= 3; j++)
      fitted_xyz[i][j] = dot(sxyz[i], R[j]) + orgi[j];

  /* rms deviation */
  temp = 0.0;
  for (i = 1; i <= n; i++) {
    ddxyz(fitted_xyz[i], exyz[i], D);
    temp += dot(D, D);
  }
  rms_value = sqrt(temp / n);

  free_dmatrix(U, 1, 3, 1, 3);
  free_dmatrix(N, 1, 4, 1, 4);
  free_dmatrix(V, 1, 4, 1, 4);

  return rms_value;
}

/* fit a plane to a set of points by least squares */
void ls_plane(double **bxyz, long n, double *pnormal, double *ppos,
              double *odist, double *adist) {
  double D[4];
  double **cov_mtx, **identityV, **V;
  long i, j, nml = 0;

  if (n < 3)
    fatal("too few atoms for least-squares fitting\n");

  cov_mtx = dmatrix(1, 3, 1, 3);
  V = dmatrix(1, 3, 1, 3);
  identityV = dmatrix(1, 3, 1, 3);

  cov_matrix(bxyz, bxyz, n, 3, cov_mtx);
  jacobi(cov_mtx, 3, D, V);

  identity_matrix(identityV, 3);
  for (i = 1; i <= 3 && !nml; i++)
    for (j = 1; j <= 3 && !nml; j++)
      if (fabs(V[i][j] - identityV[i][j]) > XEPS)
        nml = 1;
  if (nml)
    for (i = 1; i <= 3; i++)
      pnormal[i] = V[i][1];
  else { /* V is an identity matrix */
    pnormal[1] = 0.0;
    pnormal[2] = 0.0;
    pnormal[3] = 1.0;
  }

  ave_dmatrix(bxyz, n, 3, ppos);

  /* make the z-component of pnormal to be positive */
  if (pnormal[3] < 0)
    negate_xyz(pnormal);

  /* distance from the origin to the plane */
  *odist = dot(ppos, pnormal);

  /* distance from each point to the plane */
  for (i = 1; i <= n; i++)
    adist[i] = dot(bxyz[i], pnormal) - *odist;

  free_dmatrix(cov_mtx, 1, 3, 1, 3);
  free_dmatrix(V, 1, 3, 1, 3);
  free_dmatrix(identityV, 1, 3, 1, 3);
}

/* get the arbitrary rotation matrix */
void arb_rotation(double *va, double ang_deg, double **rot_mtx) {
  double c, dc, s, vlen;
  long i;

  vlen = veclen(va);
  if (vlen < XEPS) /* [0 0 0] */
    identity_matrix(rot_mtx, 3);
  else {
    for (i = 1; i <= 3; i++)
      va[i] /= vlen; /* unit vector */
    ang_deg = deg2rad(ang_deg);
    c = cos(ang_deg);
    s = sin(ang_deg);
    dc = 1 - c;
    rot_mtx[1][1] = c + dc * va[1] * va[1];
    rot_mtx[1][2] = va[1] * va[2] * dc - va[3] * s;
    rot_mtx[1][3] = va[1] * va[3] * dc + va[2] * s;
    rot_mtx[2][1] = va[1] * va[2] * dc + va[3] * s;
    rot_mtx[2][2] = c + dc * va[2] * va[2];
    rot_mtx[2][3] = va[2] * va[3] * dc - va[1] * s;
    rot_mtx[3][1] = va[1] * va[3] * dc - va[2] * s;
    rot_mtx[3][2] = va[2] * va[3] * dc + va[1] * s;
    rot_mtx[3][3] = c + dc * va[3] * va[3];
  }
}

/* angle in degrees between va and vb with vref for sign control
   va & vb are unchanged by making an additional copy of each
   all three vectors are 1-by-3 */
double vec_ang(double *va, double *vb, double *vref) {
  double ang_deg, va_cp[4], vb_cp[4];

  /* make a copy of va and vb */
  cpxyz(va, va_cp);
  cpxyz(vb, vb_cp);

  /* get orthogonal components */
  vec_orth(va_cp, vref);
  vec_orth(vb_cp, vref);

  /* angle in absolute sense */
  ang_deg = magang(va_cp, vb_cp);

  if (sign_control(va_cp, vb_cp, vref) < 0)
    ang_deg = -ang_deg;

  return ang_deg;
}

/* get the vector which has certain angle with another vector */
void get_vector(double *va, double *vref, double deg_ang, double *vo) {
  double va_cp[4], **temp;
  long i;

  cpxyz(va, va_cp); /* make a copy of va: <ana_fncs.c> */
  if (dot(va_cp, vref) > XEPS) {
    fprintf(stderr, "Angle between va/vref: %.3f degrees\n",
            magang(va_cp, vref));
    vec_orth(va_cp, vref);
  }

  temp = dmatrix(1, 3, 1, 3);

  arb_rotation(vref, deg_ang, temp);
  for (i = 1; i <= 3; i++)
    vo[i] = dot(temp[i], va_cp);
  vec_norm(vo);

  free_dmatrix(temp, 1, 3, 1, 3);
}

void rotate(double **a, long i, long j, long k, long l, double *g, double *h,
            double s, double tau) {
  *g = a[i][j];
  *h = a[k][l];
  a[i][j] = *g - s * (*h + *g * tau);
  a[k][l] = *h + s * (*g - *h * tau);
}

/* sort eigenvalues into ascending order and rearrange eigenvectors */
void eigsrt(double *d, double **v, long n) {
  double p;
  long i, j, k;

  for (i = 1; i < n; i++) {
    p = d[k = i];
    for (j = i + 1; j <= n; j++)
      if (d[j] < p)
        p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 1; j <= n; j++)
        dval_swap(&v[j][i], &v[j][k]);
    }
  }
}

void jacobi(double **a, long n, double *d, double **v) {
  long i, j, iq, ip;
  double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

  b = dvector(1, n);
  z = dvector(1, n);
  identity_matrix(v, n);
  for (ip = 1; ip <= n; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }
  for (i = 1; i <= 100; i++) {
    sm = 0.0;
    for (ip = 1; ip <= n - 1; ip++) {
      for (iq = ip + 1; iq <= n; iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm < XEPS) {
      free_dvector(z, 1, n);
      free_dvector(b, 1, n);
      eigsrt(d, v, n);
      return;
    }
    if (i < 4)
      tresh = 0.2 * sm / (n * n);
    else
      tresh = 0.0;
    for (ip = 1; ip <= n - 1; ip++) {
      for (iq = ip + 1; iq <= n; iq++) {
        g = 100.0 * fabs(a[ip][iq]);
        if (i > 4 && (fabs(d[ip]) + g) == fabs(d[ip]) &&
            (fabs(d[iq]) + g) == fabs(d[iq]))
          a[ip][iq] = 0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h = d[iq] - d[ip];
          if ((fabs(h) + g) == fabs(h))
            t = a[ip][iq] / h;
          else {
            theta = 0.5 * h / a[ip][iq];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0)
              t = -t;
          }
          c = 1.0 / sqrt(1 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq] = 0.0;
          for (j = 1; j <= ip - 1; j++)
            rotate(a, j, ip, j, iq, &g, &h, s, tau);
          for (j = ip + 1; j <= iq - 1; j++)
            rotate(a, ip, j, j, iq, &g, &h, s, tau);
          for (j = iq + 1; j <= n; j++)
            rotate(a, ip, j, iq, j, &g, &h, s, tau);
          for (j = 1; j <= n; j++)
            rotate(v, j, ip, j, iq, &g, &h, s, tau);
        }
      }
    }
    for (ip = 1; ip <= n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  fatal("too many iterations\n");
}

void dludcmp(double **a, long n, long *indx, double *d) {
  double big, dum, sum, temp;
  double *vv;
  long i, j, k;
  long imax = 0; /* initialization */

  vv = dvector(1, n);
  *d = 1.0;
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if ((temp = fabs(a[i][j])) > big)
        big = temp;
    if (big == 0.0)
      fatal("singular matrix in routine dludcmp\n");
    vv[i] = 1.0 / big;
  }
  for (j = 1; j <= n; j++) {
    for (i = 1; i < j; i++) {
      sum = a[i][j];
      for (k = 1; k < i; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i <= n; i++) {
      sum = a[i][j];
      for (k = 1; k < j; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 1; k <= n; k++)
        dval_swap(&a[imax][k], &a[j][k]);
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)
      a[j][j] = XEPS;
    if (j != n) {
      dum = 1.0 / a[j][j];
      for (i = j + 1; i <= n; i++)
        a[i][j] *= dum;
    }
  }
  free_dvector(vv, 1, n);
}

void dlubksb(double **a, long n, long *indx, double *b) {
  double sum;
  long i, ii = 0, ip, j;

  for (i = 1; i <= n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii)
      for (j = ii; j <= i - 1; j++)
        sum -= a[i][j] * b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = n; i >= 1; i--) {
    sum = b[i];
    for (j = i + 1; j <= n; j++)
      sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

void dinverse(double **a, long n, double **y) {
  double d, *col;
  long i, j, *indx;

  col = dvector(1, n);
  indx = lvector(1, n);

  dludcmp(a, n, indx, &d);

  for (j = 1; j <= n; j++) {
    for (i = 1; i <= n; i++)
      col[i] = 0.0;
    col[j] = 1.0;
    dlubksb(a, n, indx, col);
    for (i = 1; i <= n; i++)
      y[i][j] = col[i];
  }

  free_lvector(indx, 1, n);
  free_dvector(col, 1, n);
}

void get_alc_nums(char *alcname, long *num, long *nbond) {
  char str[BUF512];
  FILE *fp;

  fp = open_file(alcname, "r");

  if (fgets(str, sizeof str, fp) == NULL ||
      sscanf(str, "%ld %*s %ld", num, nbond) != 2)
    fatal("can not read atom & bond numbers: %s\n", str);

  if (!*num)
    fprintf(stderr, "ALCHEMY file <%s> has NO atoms\n", alcname);

  close_file(fp);
}

void read_alc(char *alcname, long *num, long *nbond, char **AtomName,
              double **xyz, long *ibase, long **linkage) {
  char str[BUF512], temp[BUF512];
  double c;
  long i;
  FILE *fp;

  fp = open_file(alcname, "r");

  if (fgets(str, sizeof str, fp) == NULL ||
      sscanf(str, "%ld %*s %ld", num, nbond) != 2)
    fatal("cannot read atom & bond numbers: %s\n", str);

  for (i = 1; i <= *num; i++)
    if (fgets(str, sizeof str, fp) != NULL) {
      strncpy(AtomName[i], str + 6, 2);
      AtomName[i][2] = '\0';
      upperstr(AtomName[i]); /* to upper case */
      strncpy(temp, str + 11, 9);
      temp[9] = '\0';
      if (sscanf(temp, "%lf", &xyz[i][1]) != 1)
        fatal("error reading x-coordinate\n");
      strncpy(temp, str + 20, 9);
      temp[9] = '\0';
      if (sscanf(temp, "%lf", &xyz[i][2]) != 1)
        fatal("error reading y-coordinate\n");
      strncpy(temp, str + 29, 9);
      temp[9] = '\0';
      if (sscanf(temp, "%lf", &xyz[i][3]) != 1)
        fatal("error reading z-coordinate\n");
      strncpy(temp, str + 40, 9);
      temp[9] = '\0';
      if (sscanf(temp, "%lf", &c) == 1)
        ibase[i] = lround(10.0 * c);
      else
        ibase[i] = NON_WC_IDX; /* default */
    } else
      fatal("error in reading atom records\n");

  for (i = 1; i <= *nbond; i++)
    if (fgets(str, sizeof str, fp) != NULL) {
      strncpy(temp, str + 6, 5);
      temp[5] = '\0';
      if (sscanf(temp, "%ld", &linkage[i][1]) != 1)
        fatal("error reading linkage atom 1\n");
      strncpy(temp, str + 12, 5);
      temp[5] = '\0';
      if (sscanf(temp, "%ld", &linkage[i][2]) != 1)
        fatal("error reading linkage atom 2\n");
    } else
      fatal("error in reading linkage information\n");

  close_file(fp);
}

void write_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
               long **linkage, char *alcfile) {
  char aname[5];
  long i;
  FILE *fp;

  reset_xyz(num, xyz, "f9.4");

  fp = open_file(alcfile, "w");
  fprintf(fp, "%5ld ATOMS, %5ld BONDS\n", num, nbond);
  for (i = 1; i <= num; i++) {
    if (sscanf(AtomName[i], "%s", aname) == EOF)
      strcpy(aname, "C ");
    if (isdigit((int)aname[0])) { /* like "5H1 " etc */
      aname[0] = aname[1];
      aname[1] = aname[2];
    }
    if (isdigit((int)aname[1]) || aname[1] == '\0')
      aname[1] = ' ';
    aname[2] = '\0';
    fprintf(fp, "%5ld %-2s   %9.4f%9.4f%9.4f  %9.4f\n", i, aname, xyz[i][1],
            xyz[i][2], xyz[i][3], (ibase == NULL) ? 0.0 : ibase[i] / 10.0);
  }
  for (i = 1; i <= nbond; i++)
    fprintf(fp, "%5ld %5ld %5ld\n", i, linkage[i][1], linkage[i][2]);
  close_file(fp);
}

void free_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
              long zero_1, long **linkage) {
  free_cmatrix(AtomName, 1, num, 0, 2);
  free_dmatrix(xyz, 1, num, 1, 3);
  free_lvector(ibase, 1, num);
  free_lmatrix(linkage, 1, nbond, zero_1, 2);
}

/* sort a long vector into ascending order & keep the index
   use Shell's method as in NR in C book Ed. 2, pp.331-332 */
void lsort(long n, long *a, long *idx) {
  long v;
  long i, inc, iv, j;

  inc = 1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);

  for (i = 1; i <= n; i++)
    idx[i] = i;

  do {
    inc /= 3;
    for (i = inc + 1; i <= n; i++) {
      v = a[i];
      iv = idx[i];
      j = i;
      while (a[j - inc] > v) {
        a[j] = a[j - inc];
        idx[j] = idx[j - inc];
        j -= inc;
        if (j <= inc)
          break;
      }
      a[j] = v;
      idx[j] = iv;
    }
  } while (inc > 1);
}

/* reverse a long vector: output replaces the original */
void lreverse(long ia, long n, long *lvec) {
  long i, *ltmp;

  ltmp = lvector(1, n);
  for (i = 1; i <= n; i++)
    ltmp[i] = lvec[n + ia - i];
  for (i = 1; i <= n; i++)
    lvec[ia + i - 1] = ltmp[i];

  free_lvector(ltmp, 1, n);
}

/* write a title to XFIG 3.2 format */
void fig_title(FILE *fp) {
  fprintf(fp,
          "#FIG 3.2  Creator: %s\n"
          "Portrait\n"
          "Flush Left\n"
          "Inches\n"
          "Letter\n"
          "100.00\n"
          "Single\n"
          "-2\n"
          "1200 2\n\n",
          Gvars.X3DNA_VER);
}

void ps_title_cmds(FILE *fp, char *imgfile, long *bbox) {
  char BDIR[BUF1K], str[BUF512];
  char *ps_image_par = "ps_image.par";
  long i;
  time_t run_time;
  FILE *fpp;

  run_time = time(NULL);

  fprintf(fp, "%%!PS-Adobe-3.0\n");
  fprintf(fp, "%%%%Title: (%s)\n", imgfile);
  fprintf(fp, "%%%%Creator: (%s)\n", Gvars.X3DNA_VER);
  strcpy(str, ctime(&run_time));
  str[strlen(str) - 1] = '\0'; /* change \n to \0 */
  fprintf(fp, "%%%%CreationDate: (%s)\n", str);
  fprintf(fp, "%%%%Orientation: Portrait\n");
  fprintf(fp, "%%%%BoundingBox: ");
  for (i = 1; i <= 4; i++)
    fprintf(fp, "%6ld", bbox[i]);
  fprintf(fp, "\n\n");

  /* read in color parameter file */
  get_BDIR(BDIR, ps_image_par);
  strcat(BDIR, ps_image_par);
  fprintf(stderr, " ...... reading file: %s ...... \n", ps_image_par);

  fpp = open_file(BDIR, "r");
  while (fgets(str, sizeof str, fpp) != NULL)
    fprintf(fp, "%s", str);
  close_file(fpp);
}

/* reset x/y coordinates to XFIG units */
void get_fig_xy(long num, double **xyz, long nO, double **oxyz, long *urxy,
                long frame_box, FILE *fp) {
  char *format = "%6ld%6ld";
  double paper_size[2] = {8.5, 11.0}; /* US letter size */
  double max_xy[3];
  long i, j, llxy[3];

  /* change y-axis according to XFIG coordinate system */
  for (i = 1; i <= num; i++)
    xyz[i][2] = urxy[2] - xyz[i][2];
  if (nO)
    for (i = 1; i <= nO; i++)
      oxyz[i][2] = urxy[2] - oxyz[i][2];

  max_dmatrix(xyz, num, 2, max_xy);
  for (i = 1; i <= 2; i++)
    urxy[i] = lround(max_xy[i]); /* new urxy */

  /* 1/4 offset the figure on a US letter (8.5in-by-11in) */
  for (i = 1; i <= 2; i++)
    llxy[i] = lround(0.25 * (paper_size[i - 1] * 1200 - urxy[i]));
  for (i = 1; i <= num; i++)
    for (j = 1; j <= 2; j++)
      xyz[i][j] += llxy[j];
  if (nO)
    for (i = 1; i <= nO; i++)
      for (j = 1; j <= 2; j++)
        oxyz[i][j] += llxy[j];

  max_dmatrix(xyz, num, 2, max_xy);
  for (i = 1; i <= 2; i++)
    urxy[i] = lround(max_xy[i]);

  fig_title(fp);

  if (frame_box) {
    /* draw a box around the figure at the deepest depth (999) */
    fprintf(fp, "# draw a boundary box here\n");
    fprintf(fp, "2 3 0 1 0 0 999 0 -1 0.0 2 1 0 0 0 5\n");
    fprintf(fp, format, llxy[1] - FIG_BOUND, llxy[2] - FIG_BOUND);
    fprintf(fp, format, urxy[1] + FIG_BOUND, llxy[2] - FIG_BOUND);
    fprintf(fp, format, urxy[1] + FIG_BOUND, urxy[2] + FIG_BOUND);
    fprintf(fp, format, llxy[1] - FIG_BOUND, urxy[2] + FIG_BOUND);
    fprintf(fp, format, llxy[1] - FIG_BOUND, llxy[2] - FIG_BOUND);
    fprintf(fp, "\n\n");
  }
}

/* get the point position of x-, y-, or z-axis */
void get_pxy(double *xy1, double *xy2, double r, double *px, double *py) {
  double dx, dy, dd;

  dx = xy2[1] - xy1[1];
  dy = xy2[2] - xy1[2];
  dd = sqrt(dx * dx + dy * dy);

  if (dd > XEPS) {
    *px = xy2[1] + r * dx / dd;
    *py = xy2[2] + r * dy / dd;
  } else {
    *px = xy2[1];
    *py = xy2[2];
  }
}

void alc2fig(long nobj, long *idx, long *depth, long **allobj, double **blkxyz,
             double **oxyz, long *ibase, long faces[][5], long *opts,
             FILE *fp) {
  char *format = "%6.0f%6.0f";
  double dot_sep, msat, Msat, px, py;
  long is_color, same_faces, updown, mgroove;
  long dlcol, dwidth, w1, w2, line_width, join_style, cap_style, o_sides, mfcol;
  long bcol_code, i, ib, ie, ioffset8, ip, j, k, Mfill, mfill;
  long **bc_idx; /* [2][7] -- color vs black/white, with 7 types total */

  is_color = opts[2]; /* decomposed for clarity */
  same_faces = opts[7];
  updown = opts[9];
  mgroove = opts[10];

  bc_idx = lmatrix(0, 1, 0, 6); /* [2][7] */
  get_fig_pars(&dot_sep, &dlcol, &dwidth, &w1, &w2, bc_idx, &msat, &Msat,
               &o_sides, &line_width, &join_style, &cap_style, &mfcol);

  /* render each object from inside to outside */
  for (i = 1; i <= nobj; i++) {
    j = idx[i];             /* object number */
    if (allobj[1][j] > 0) { /* blocks */
      bcol_code = bc_idx[is_color][ibase[allobj[1][j]]];
      k = lround(msat * 20);
      if (mgroove) /* could still be in color */
        mfill = 20 - k;
      else
        mfill = (bcol_code) ? 20 + k : 20 - k;
      k = lround(Msat * 20);
      if (mgroove) /* could still be in color */
        Mfill = k;
      else
        Mfill = (bcol_code) ? 40 - k : k;

      if (!same_faces) {
        if ((!updown && !allobj[2][j]) ||  /* minor groove side */
            (updown && allobj[2][j] == 4)) /* upper face */
          fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 %2ld", line_width, bcol_code,
                  bcol_code, depth[i], mfill);
        else if ((!updown && allobj[2][j] == 1) || /* major groove side: back */
                 (updown && allobj[2][j] == 5))    /* lower face */
          fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 %2ld", line_width, bcol_code,
                  bcol_code, depth[i], Mfill);
        else {           /* other sides */
          if (bcol_code) /* color */
            fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 20", line_width, bcol_code,
                    o_sides, depth[i]);
          else /* black & white */
            fprintf(fp, "2 3 0 %2ld 0 0 %4ld 0 0", line_width, depth[i]);
        }
      } else {
        if (mgroove && !allobj[2][j]) /* minor groove side */
          fprintf(fp, "2 3 0 %2ld 0 %2ld %4ld 0 %2ld", line_width, bc_idx[0][0],
                  depth[i], mfcol);
        else
          fprintf(fp, "2 3 0 %2ld 0 %2ld %4ld 0 %2ld", line_width, bcol_code,
                  depth[i], (bcol_code) ? mfill : Mfill);
      }

      fprintf(fp, " 0.0 %2ld %2ld 0 0 0 5\n", join_style, cap_style);
      ioffset8 = (allobj[1][j] - 1) * 8;
      for (k = 0; k < 5; k++) {
        ip = ioffset8 + faces[allobj[2][j]][k];
        fprintf(fp, format, blkxyz[ip][1], blkxyz[ip][2]);
      }
      fprintf(fp, "\n");
    } else { /* O connection lines */
      if (!is_color)
        dlcol = 0;         /* black */
      if (!allobj[1][j]) { /* origin lines */
        fprintf(fp, "2 1 2 %2ld %2ld 0 %4ld 0 -1 %5.1f", dwidth, dlcol,
                depth[i], dot_sep);
        fprintf(fp, " %2ld %2ld 0 0 0 2\n", join_style, cap_style);
      } else if (allobj[1][j] == -1) { /* helix axis */
        fprintf(fp, "2 1 0 %2ld %2ld 0 %4ld 0 -1 0.0", w2, dlcol, depth[i]);
        fprintf(fp, " %2ld %2ld 0 0 0 2\n", join_style, cap_style);
      } else { /* reference axis */
        fprintf(fp, "2 1 0 %2ld %2ld 0 %4ld 0 -1 0.0", w1, dlcol, depth[i]);
        fprintf(fp, " %2ld %2ld 0 1 0 2\n", join_style, cap_style);
        fprintf(fp, "  2 0 2 36 66\n");
      }
      ib = allobj[2][j] / 10000;
      fprintf(fp, format, oxyz[ib][1], oxyz[ib][2]);
      ie = allobj[2][j] % 10000;
      fprintf(fp, format, oxyz[ie][1], oxyz[ie][2]);
      fprintf(fp, "\n");
      if (allobj[1][j] <= -2) {
        get_pxy(oxyz[ib], oxyz[ie], 10.0, &px, &py);
        fprintf(fp, "4 1 %2ld %4ld 0 18 18 0.0 4 165 165", dlcol, depth[i]);
        /* Helvetica-Bold: #18 */
        fprintf(fp, format, px, py);
        if (allobj[1][j] == -2)
          fprintf(fp, " x\\001\n");
        else if (allobj[1][j] == -3)
          fprintf(fp, " y\\001\n");
        else
          fprintf(fp, " z\\001\n");
      }
    }
  }

  free_lmatrix(bc_idx, 0, 1, 0, 6);
}

/* reset x/y coordinates to PS units */
void get_ps_xy(char *imgfile, long *urxy, long frame_box, FILE *fp) {
  char *format = "%6ld%6ld";
  double paper_size[2] = {8.5, 11.0}; /* US letter size */
  long i;
  long bbox[5], llxy[3];

  /* centralize the figure on a US letter (8.5in-by-11in) */
  for (i = 1; i <= 2; i++)
    llxy[i] = lround(0.5 * (paper_size[i - 1] * 72 - urxy[i]));

  /* boundary box */
  for (i = 1; i <= 2; i++) {
    bbox[i] = llxy[i] - PS_BOUND;
    bbox[i + 2] = urxy[i] + llxy[i] + PS_BOUND;
  }

  ps_title_cmds(fp, imgfile, bbox);

  fprintf(fp, "%6ld%6ld translate\n\n", llxy[1], llxy[2]);

  if (frame_box) {
    /* draw a box around the figure */
    fprintf(fp, "NP ");
    fprintf(fp, format, -PS_BOUND, -PS_BOUND);
    fprintf(fp, format, urxy[1] + PS_BOUND, -PS_BOUND);
    fprintf(fp, format, urxy[1] + PS_BOUND, urxy[2] + PS_BOUND);
    fprintf(fp, format, -PS_BOUND, urxy[2] + PS_BOUND);
    fprintf(fp, " DB stroke\n\n");
  }
}

void alc2ps(long nobj, long *idx, long **allobj, double **blkxyz, double **oxyz,
            long *ibase, long faces[][5], long *opts, FILE *fp) {
  char bname;
  char *format = "%7.1f%7.1f";
  char bc_idx[2][8] = {
      "XXXXXXX", /* black & white style: [8] to allow for '\0' */
      CX_LIST};  /* color-coded for ACGITUX */
  double px, py;
  long is_color, same_faces, updown, mgroove;
  long i, ib, ie, ioffset8, ip, j, k;

  is_color = opts[2]; /* decomposed for clarity */
  same_faces = opts[7];
  updown = opts[9];
  mgroove = opts[10];

  /* render each object from inside to outside */
  for (i = 1; i <= nobj; i++) {
    j = idx[i]; /* object number */
    fprintf(fp, "NP ");
    if (allobj[1][j] > 0) { /* blocks */
      bname = bc_idx[is_color][ibase[allobj[1][j]]];
      fprintf(fp, "%cl ", (!same_faces) ? bname : 'X');
      ioffset8 = (allobj[1][j] - 1) * 8;
      for (k = 0; k < 4; k++) {
        ip = ioffset8 + faces[allobj[2][j]][k];
        fprintf(fp, format, blkxyz[ip][1], blkxyz[ip][2]);
      }
      fprintf(fp, " DB\n");
      if (!same_faces) {
        if ((!updown && !allobj[2][j]) ||  /* minor groove side */
            (updown && allobj[2][j] == 4)) /* upper face */
          fprintf(fp, "  gsave %cm grestore stroke\n", bname);
        else if ((!updown && allobj[2][j] == 1) || /* major groove side: back */
                 (updown && allobj[2][j] == 5))    /* lower face */
          fprintf(fp, "  gsave %cM grestore stroke\n", bname);
        else {              /* other sides */
          if (bname == 'X') /* black & white */
            fprintf(fp, "  gsave 1.0 FB grestore stroke\n");
          else
            fprintf(fp, "  gsave OTHER_SIDES grestore stroke\n");
        }
      } else {
        if (mgroove && !allobj[2][j]) /* minor groove side */
          fprintf(fp, "  gsave Sm grestore stroke\n");
        else {
          if (bname == 'X') /* black & white */
            fprintf(fp, "  gsave XM grestore stroke\n");
          else
            fprintf(fp, "  gsave %cm grestore stroke\n", bname);
        }
      }
    } else { /* O connection lines */
      (is_color) ? fprintf(fp, "Dl ") : fprintf(fp, "Xl ");
      ib = allobj[2][j] / 10000;
      fprintf(fp, format, oxyz[ib][1], oxyz[ib][2]);
      ie = allobj[2][j] % 10000;
      fprintf(fp, format, oxyz[ie][1], oxyz[ie][2]);
      if (!allobj[1][j]) /* origin lines */
        fprintf(fp, "  gsave Dw Ds LN grestore\n");
      else if (allobj[1][j] == -1) /* helix axis */
        fprintf(fp, "  gsave W2 LN grestore\n");
      else { /* reference axis */
        fprintf(fp, "  gsave W1 LN grestore\n");

        get_pxy(oxyz[ib], oxyz[ie], 8.0, &px, &py);
        fprintf(fp, format, px, py);
        fprintf(fp, "  moveto");
        if (allobj[1][j] == -2)
          fprintf(fp, " (x) SCENTER\n");
        else if (allobj[1][j] == -3)
          fprintf(fp, " (y) SCENTER\n");
        else
          fprintf(fp, " (z) SCENTER\n");
      }
    }
  }
  fprintf(fp, "\nshowpage\n");
}

/* get base ring atom index in one residue */
void bring_atoms(long ib, long ie, long ra_num, char **AtomName, long *nmatch,
                 long *batom) {
  static char *RingAtom[] = {RA_LIST};
  long i, j;

  *nmatch = 0;

  for (i = 0; i < ra_num; i++) {
    j = find_1st_atom(RingAtom[i], AtomName, ib, ie, "in base ring atoms");
    if (j)
      batom[++*nmatch] = j;
  }
}

/* get base ring atom index for all residues: num_ring */
void all_bring_atoms(long num_residue, long *RY, long **seidx, char **AtomName,
                     long *num_ring, long **ring_atom) {
  long i, j, nmatch;

  *num_ring = 0;
  for (i = 1; i <= num_residue; i++) {
    if (RY[i] < 0) { /* non-base residue */
      ring_atom[i][10] = -1;
      continue;
    }
    j = (RY[i] == 1) ? 9 : 6;
    bring_atoms(seidx[i][1], seidx[i][2], j, AtomName, &nmatch, ring_atom[i]);
    if (nmatch == j) {
      ring_atom[i][10] = j;
      ++*num_ring;
    }
  }
}

/* get base index for coloring purpose */
void base_idx(long num, char *bseq, long *ibase, long single) {
  static char *cmn_base = CB_LIST, *pchar;
  long i;

  if (single) { /* for a single case */
    if ((pchar = strchr(cmn_base, toupper((int)*bseq))) != NULL)
      *ibase = pchar - cmn_base;
    else
      *ibase = NON_WC_IDX;
  } else {
    for (i = 1; i <= num; i++)
      if ((pchar = strchr(cmn_base, toupper((int)bseq[i]))) != NULL)
        ibase[i] = pchar - cmn_base;
      else
        ibase[i] = NON_WC_IDX;
  }
}

/* get basepair index for coloring purpose */
long basepair_idx(char *bpi) {
  static char *WC[9] = {WC_LIST};
  long i, bidx;

  i = find_1st_atom(bpi, WC, 1, 8, "");
  if (!i)
    bidx = NON_WC_IDX; /* non A-T/G-C pair */
  else if (i <= 2)     /* AT or AU */
    bidx = 0;          /* 0 for A */
  else if (i <= 4)     /* TA or UA */
    bidx = 4;          /* 4 for T */
  else if (i <= 6)     /* GC or IC */
    bidx = 2;          /* 2 for G */
  else                 /* CG or CI */
    bidx = 1;          /* 1 for C */

  return bidx;
}

/* given plane normal and its center, project all coordinates onto it */
void plane_xyz(long num, double **xyz, double *ppos, double *nml,
               double **nxyz) {
  long i, j;
  double temp, d[4];

  for (i = 1; i <= num; i++) {
    ddxyz(ppos, xyz[i], d);
    temp = dot(d, nml);
    for (j = 1; j <= 3; j++)
      nxyz[i][j] = ppos[j] + d[j] - temp * nml[j];
  }
}

/* project base atoms onto its least-squares plane defined by ring atoms */
void prj2plane(long num, long ra_num, char **AtomName, double **xyz, double z0,
               double **nxyz) {
  double ang, temp;
  double zaxis[4] = {EMPTY_NUMBER, 0.0, 0.0, 1.0};
  double adist[10], hinge[4], ppos[4], z[4];
  double **bxyz, **rmtx;

  long i, j, nmatch;
  long batom[10];

  /* find base ring atoms */
  bring_atoms(1, num, ra_num, AtomName, &nmatch, batom);

  /* base least-squares plane */
  bxyz = dmatrix(1, nmatch, 1, 3);

  for (i = 1; i <= nmatch; i++)
    cpxyz(xyz[batom[i]], bxyz[i]);
  ls_plane(bxyz, nmatch, z, ppos, &temp, adist);

  /* get the new set of coordinates */
  plane_xyz(num, xyz, ppos, z, nxyz);

  /* reorient the structure to make z-coordinate zero */
  rmtx = dmatrix(1, 3, 1, 3);
  if (z0) {
    cross(z, zaxis, hinge);
    ang = magang(z, zaxis);
    arb_rotation(hinge, ang, rmtx);
    for (i = 1; i <= num; i++)
      for (j = 1; j <= 3; j++)
        nxyz[i][j] = dot(xyz[i], rmtx[j]);
    for (i = 1; i <= num; i++)
      nxyz[i][3] -= nxyz[1][3];
  }
  free_dmatrix(bxyz, 1, nmatch, 1, 3);
  free_dmatrix(rmtx, 1, 3, 1, 3);
}

/* reset x & y coordinates to fit the scale */
void adjust_xy(long num, double **xyz, long nO, double **oxyz,
               double scale_factor, long default_size, long *urxy) {
  long i, j;
  double temp;
  double max_xy[3], min_xy[3];

  max_dmatrix(xyz, num, 2, max_xy);
  min_dmatrix(xyz, num, 2, min_xy);

  /* get maximum dx or dy */
  temp = dval_max(max_xy[1] - min_xy[1], max_xy[2] - min_xy[2]);

  scale_factor = fabs(scale_factor);
  if (scale_factor < XEPS)
    scale_factor = default_size / temp;

  fprintf(stderr, "\n ...... scale factor: %.2f ...... \n", scale_factor);

  move_position(xyz, num, 2, min_xy);
  for (i = 1; i <= num; i++)
    for (j = 1; j <= 2; j++)
      xyz[i][j] *= scale_factor;
  if (nO) {
    move_position(oxyz, nO, 2, min_xy);
    for (i = 1; i <= nO; i++)
      for (j = 1; j <= 2; j++)
        oxyz[i][j] *= scale_factor;
  }
  max_dmatrix(xyz, num, 2, max_xy);
  for (i = 1; i <= 2; i++)
    urxy[i] = lround(max_xy[i]);
}

void get_depth(long nobj, long *zval, long *depth) {
  double temp;
  long depth_low = 991, depth_up = 11; /* depth level */
  long i, j;

  /* reset zval to [depth_low -- depth_up] for depth level */
  j = depth_low - depth_up;
  temp = zval[nobj] - zval[1];
  for (i = 1; i <= nobj; i++)
    depth[i] = lround(depth_low - j * (zval[i] - zval[1]) / temp);
}

void raster3d_header(long num, double **xyz, double scale_factor,
                     long no_header, long frame_box, FILE *fp) {
  char BDIR[BUF1K], str[BUF512], *header_file = "my_header.r3d";
  double temp, ave_xyz[4], min_xyz[4], max_xyz[4];
  double rad = 0.06, rgbv[4] = {0.0, 0.25, 0.25, 0.25};
  long i, itype = 3;
  FILE *fpp;

  min_dmatrix(xyz, num, 3, min_xyz);
  max_dmatrix(xyz, num, 3, max_xyz);
  avexyz(max_xyz, min_xyz, ave_xyz);
  temp = dval_max(max_xyz[1] - min_xyz[1], max_xyz[2] - min_xyz[2]);

  if (!no_header) { /* write header section */
    if (scale_factor < XEPS)
      scale_factor = 6.0 + temp;
    fprintf(stderr, "\n ...... scale factor: %.2f ...... \n", scale_factor);

    get_BDIR(BDIR, header_file);
    strcat(BDIR, header_file);
    fprintf(stderr, " ...... reading file: %s ...... \n", header_file);

    fpp = open_file(BDIR, "r");
    for (i = 1; i <= 20; i++) {
      if (fgets(str, sizeof str, fpp) == NULL)
        fatal("error reading header.r3d\n");
      (i != 16) ? fputs(str, fp)
                : fprintf(fp, "%9.3f%9.3f%9.3f%9.3f\n", -ave_xyz[1],
                          -ave_xyz[2], -ave_xyz[3], scale_factor);
    }
    close_file(fpp);
  }

  if (frame_box) {
    fprintf(fp, "### Section of the frame box: 4 lines\n");
    ave_xyz[3] = min_xyz[3] = max_xyz[3];
    /* box strongly influenced by z-coordinate */
    ave_xyz[1] = max_xyz[1]; /* lower-right corner */
    ave_xyz[2] = min_xyz[2];
    r3d_rod(itype, min_xyz, ave_xyz, rad, rgbv, fp);
    r3d_rod(itype, ave_xyz, max_xyz, rad, rgbv, fp);
    ave_xyz[1] = min_xyz[1]; /* upper-left corner */
    ave_xyz[2] = max_xyz[2];
    r3d_rod(itype, min_xyz, ave_xyz, rad, rgbv, fp);
    r3d_rod(itype, ave_xyz, max_xyz, rad, rgbv, fp);
  }
}

/* read in parameters for Raster3D input */
void get_r3dpars(double **base_col, double *hb_col, double *width3,
                 double **atom_col, char *label_style) {
  char BDIR[BUF1K], str[BUF512], *raster3d_par = "raster3d.par";
  char *format = "%lf %lf %lf", *format4 = "%lf %lf %lf %lf";
  long i;
  FILE *fp;

  get_BDIR(BDIR, raster3d_par);
  strcat(BDIR, raster3d_par);
  fp = open_file(BDIR, "r");
  if (Gvars.VERBOSE)
    fprintf(stderr, " ...... reading file: %s ...... \n", raster3d_par);

  if (fgets(str, sizeof str, fp) == NULL) /* skip one line */
    fatal("error in reading comment line\n");
  for (i = 0; i <= NBASECOL; i++)
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format, &base_col[i][1], &base_col[i][2],
               &base_col[i][3]) != 3)
      fatal("error reading base residue RGB color\n");
  if (fgets(str, sizeof str, fp) == NULL) /* skip one line */
    fatal("error in reading comment line\n");
  if (fgets(str, sizeof str, fp) == NULL ||
      sscanf(str, format4, &hb_col[1], &hb_col[2], &hb_col[3], &hb_col[4]) != 4)
    fatal("error reading H-bond RGB color\n");
  if (fgets(str, sizeof str, fp) == NULL) /* skip one line */
    fatal("error in reading comment line\n");
  if (fgets(str, sizeof str, fp) == NULL ||
      sscanf(str, "%lf %lf %lf", &width3[1], &width3[2], &width3[3]) != 3)
    fatal("error cylinder radius for bp-center line & bp 1 & 2\n");
  if (fgets(str, sizeof str, fp) == NULL) /* skip one line */
    fatal("error in reading comment line\n");
  for (i = 0; i <= NATOMCOL; i++)
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format, &atom_col[i][1], &atom_col[i][2],
               &atom_col[i][3]) != 3)
      fatal("error reading atom RGB color\n");
  if (fgets(str, sizeof str, fp) == NULL) /* skip one line */
    fatal("error in reading comment line\n");
  if (fgets(str, sizeof str, fp) == NULL)
    fatal("error reading label style\n");
  strcpy(label_style, str);

  close_file(fp);
}

/* write a record of round-ended cylinder (itype = 3) or flat-ended (5) or
 * comments */
void r3d_rod(long itype, double *xyz1, double *xyz2, double rad, double *rgbv,
             FILE *fp) {
  static char *format = "%9.3f";
  long i;

  if (itype == 3 || itype == 5)
    fprintf(fp, "%ld\n", itype);
  else
    fprintf(fp, "#5\n#"); /* comments and default to type 5 */
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, xyz1[i]);
  fprintf(fp, format, rad);
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, xyz2[i]);
  fprintf(fp, format, rad);
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, rgbv[i]);
  fprintf(fp, "\n");
}

/* draw dashed line for H-bonds, and lines connecting base-pair centers */
void r3d_dash(double *xyz1, double *xyz2, double hb_width, double *hb_col,
              FILE *fp) {
  double distance, dnum, dxyz[4], m1[4], m2[4];
  long i, j, num;

  ddxyz(xyz1, xyz2, dxyz);
  distance = veclen(dxyz);
  vec_norm(dxyz); /* normal vector */

  num = lround(hb_col[4] * distance);
  if (num % 2 == 0)
    num++;      /* odd "num" for ending on last point */
  if (num <= 1) /* solid line */
    r3d_rod(5, xyz1, xyz2, hb_width, hb_col, fp);
  else {
    dnum = (double)num;
    for (i = 0; i < num; i += 2) {
      for (j = 1; j <= 3; j++) {
        m1[j] = xyz1[j] + dxyz[j] * distance * i / dnum;
        m2[j] = xyz1[j] + dxyz[j] * distance * (i + 1) / dnum;
      }
      r3d_rod(5, m1, m2, hb_width, hb_col, fp);
    }
  }
}

/* write a record of sphere */
void r3d_sphere(double *xyz1, double rad, double *rgbv, FILE *fp) {
  static char *format = "%9.3f";
  long i;

  fprintf(fp, "2\n");
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, xyz1[i]);
  fprintf(fp, format, rad);
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, rgbv[i]);
  fprintf(fp, "\n");
}

void cpk_model(long num, long *idx, double **xyz, double ballrad,
               double **colrgb, FILE *fp) {
  double vdw_radii[NELE];
  long i;

  if (ballrad <= 0.0)
    return;

  atom_info(3, NULL, NULL, vdw_radii);
  fprintf(fp, "###\n### The following section is for ball/CPK model\n");
  for (i = 1; i <= num; i++)
    r3d_sphere(xyz[i], vdw_radii[idx[i]] * ballrad, colrgb[i], fp);
}

/* write a record of triangle (itype = 1, default) or plane (itype = 6) */
void r3d_tripln(long itype, double *xyz1, double *xyz2, double *xyz3,
                double *rgbv, FILE *fp) {
  static char *format = "%9.3f";
  long i;

  fprintf(fp, "%d\n", (itype != 6) ? 1 : 6);
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, xyz1[i]);
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, xyz2[i]);
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, xyz3[i]);
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, rgbv[i]);
  fprintf(fp, "\n");
}

void r3d_block_edge(double *rgbv, long ioffset8, double **blkxyz, double w1,
                    FILE *fp) {
  static long blk_lkg[12][2] = {{1, 2}, {2, 3}, {3, 4}, {1, 4}, {1, 5}, {2, 6},
                                {3, 7}, {4, 8}, {5, 6}, {6, 7}, {7, 8}, {5, 8}};
  long j;

  for (j = 0; j < 12; j++)
    r3d_rod(3, blkxyz[ioffset8 + blk_lkg[j][0]],
            blkxyz[ioffset8 + blk_lkg[j][1]], w1, rgbv, fp);
}

/* label the base in the center of six-membered ring */
void base_label(double **rxyz, char *label_style, double *rgbv, char *bname_num,
                FILE *fp) {
  static char *format = "%9.3f";
  double cxyz[4];
  long i;

  avexyz(rxyz[1], rxyz[4], cxyz);     /* N1 + C4 */
  fprintf(fp, "10\n%s", label_style); /* label_style has \n */
  fprintf(fp, "11\n");
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, cxyz[i]);
  for (i = 1; i <= 3; i++)
    fprintf(fp, format, rgbv[i]);
  fprintf(fp, "\n%s\n", bname_num);
}

void fill_base_ring(long num_residue, long num_ring, long **ring_atom,
                    double **xyz, long *ibase, char *bseq, double **base_col,
                    char *label_style, long label_ring, long *ResSeq,
                    FILE *fp) {
  char bname_num[BUF512];
  double **rxyz;
  long i, j;

  rxyz = dmatrix(1, 9, 1, 3); /* base ring atoms */

  fprintf(fp, "###\n### The following section is for %ld filled base rings\n",
          num_ring);
  for (i = 1; i <= num_residue; i++) {
    if (ring_atom[i][10] <= 0)
      continue;
    for (j = 1; j <= ring_atom[i][10]; j++)
      cpxyz(xyz[ring_atom[i][j]], rxyz[j]);
    r3d_tripln(1, rxyz[3], rxyz[4], rxyz[6], base_col[ibase[i]], fp);
    r3d_tripln(1, rxyz[1], rxyz[2], rxyz[3], base_col[ibase[i]], fp);
    r3d_tripln(1, rxyz[4], rxyz[5], rxyz[6], base_col[ibase[i]], fp);
    r3d_tripln(1, rxyz[1], rxyz[3], rxyz[6], base_col[ibase[i]], fp);
    if (ring_atom[i][10] == 9) {
      r3d_tripln(1, rxyz[6], rxyz[7], rxyz[8], base_col[ibase[i]], fp);
      r3d_tripln(1, rxyz[1], rxyz[8], rxyz[9], base_col[ibase[i]], fp);
      r3d_tripln(1, rxyz[1], rxyz[6], rxyz[8], base_col[ibase[i]], fp);
    }
    if (label_ring) {
      sprintf(bname_num, "%c%ld", bseq[i], ResSeq[ring_atom[i][1]]);
      base_label(rxyz, label_style, base_col[ibase[i]], bname_num, fp);
    }
  }

  free_dmatrix(rxyz, 1, 9, 1, 3);
}

void process_alc(char *alcfile, char *imgfile, double scale_factor,
                 long *opts) {
  char **AtomName;
  double **blkxyz, **oxyz, **xyz;
  long i, j, k, nbond, nobj, npoint_blk, num, num_blk = 0, num_blk8;
  long ioffset8, ioffset_blk, nN = 0, nO = 0, nO_all = 0, nO_lkg = 0, nO_tmp;
  long Nidx[6], *blkibase, *ibase, **linkage;

  /* read in ALCHEMY file */
  get_alc_nums(alcfile, &num, &nbond);
  AtomName = cmatrix(1, num, 0, 2);
  xyz = dmatrix(1, num, 1, 3);
  ibase = lvector(1, num);
  linkage = lmatrix(1, nbond, 0, 2);
  read_alc(alcfile, &num, &nbond, AtomName, xyz, ibase, linkage);

  /* assume the first 8 atoms are 4 Ns followed by 4 Cs in
     a specific numbering scheme */
  for (i = 1; i <= num; i++)
    if (AtomName[i][0] == 'N')
      ++nN;
  if (!nN || nN % 4)
    fatal("wrong type of standard block file\n");
  else
    num_blk = nN / 4;

  nN = 0;
  for (i = 1; i <= num && nN < 5; i++)
    if (AtomName[i][0] == 'N')
      Nidx[++nN] = i;
  npoint_blk = (num_blk > 1) ? Nidx[5] - Nidx[1] : num;

  /* get the xyz coordinates of 8-point-per-block for each */
  num_blk8 = num_blk * 8;
  blkxyz = dmatrix(1, num_blk8, 1, 3);
  blkibase = lvector(1, num_blk);
  for (i = 1; i <= num_blk; i++) {
    ioffset8 = (i - 1) * 8;
    ioffset_blk = (i - 1) * npoint_blk;
    blkibase[i] = ibase[ioffset_blk + 1];
    if (!lval_in_range(blkibase[i], 0, 5)) {
      fprintf(stderr, "uncommon base. set to default color coding\n");
      blkibase[i] = NON_WC_IDX;
    }
    for (j = 1; j <= 8; j++)
      cpxyz(xyz[ioffset_blk + j], blkxyz[ioffset8 + j]);
  }
  nobj = num_blk * 6; /* blocks */

  /* get the total number of O atoms:
     "O ", "OH", "OO", "OX", "OY", "OZ" */
  for (i = 1; i <= num; i++)
    if (AtomName[i][0] == 'O')
      nO_all++;
  oxyz = (nO_all) ? dmatrix(1, nO_all, 1, 3) : NULL;

  if (opts[3]) { /* draw a line to link neighbor origins */
    nO_tmp = 0;
    for (i = 1; i <= num; i++)
      if (!strcmp(AtomName[i], "O ")) {
        ++nO_tmp;
        cpxyz(xyz[i], oxyz[nO_tmp]);
      }
    if (nO_tmp) {
      k = nO_tmp - 1;
      for (i = 1; i <= k; i++) {
        linkage[i][0] = 0;
        linkage[i][1] = i;
        linkage[i][2] = i + 1;
      }
      nO_lkg = k;
      nO = nO_tmp;
    }
  }
  if (opts[5]) { /* draw helix-axis */
    nO_tmp = 0;
    for (i = 1; i <= num; i++)
      if (!strcmp(AtomName[i], "OH")) {
        ++nO_tmp;
        cpxyz(xyz[i], oxyz[nO + nO_tmp]);
      }
    if (nO_tmp) {
      if (nO_tmp % 2)
        fatal("wrong number of helix axis origins\n");
      k = nO_tmp / 2;
      for (i = 1; i <= k; i++) {
        j = nO_lkg + i;
        linkage[j][0] = -1;
        linkage[j][1] = nO + 2 * i - 1;
        linkage[j][2] = nO + 2 * i;
      }
      nO_lkg += k;
      nO += nO_tmp;
    }
  }
  if (opts[6]) { /* draw global reference axis */
    nO_tmp = 0;
    for (i = 1; i <= num; i++)
      if (!strcmp(AtomName[i], "OO") || !strcmp(AtomName[i], "OX") ||
          !strcmp(AtomName[i], "OY") || !strcmp(AtomName[i], "OZ")) {
        ++nO_tmp;
        cpxyz(xyz[i], oxyz[nO + nO_tmp]);
      }
    if (nO_tmp) {
      k = nO_tmp - 1;
      for (i = 1; i <= k; i++) {
        j = nO_lkg + i;
        linkage[j][0] = -(i + 1);
        linkage[j][1] = nO + 1;
        linkage[j][2] = nO + i + 1;
      }
      nO_lkg += k;
      nO += nO_tmp;
    }
  }

  alc_3images(opts, nobj, num_blk, num_blk8, nO_lkg, nO, oxyz, blkxyz, blkibase,
              linkage, scale_factor, imgfile);

  if (oxyz != NULL)
    free_dmatrix(oxyz, 1, nO_all, 1, 3);
  free_alc(num, nbond, AtomName, xyz, ibase, 0, linkage);
  free_dmatrix(blkxyz, 1, num_blk8, 1, 3);
  free_lvector(blkibase, 1, num_blk);
}

/* From ALCHEMY format to three image formats: PS, XFIG, Raster3D */
void alc_3images(long *opts, long nobj, long num_blk, long num_blk8,
                 long nO_lkg, long nO, double **oxyz, double **blkxyz,
                 long *blkibase, long **linkage, double scale_factor,
                 char *imgfile) {
  char label_style[BUF512];
  double width3[4], hb_col[5], *pd, **atom_col, **base_col;
  static long faces[6][5] = {
      {1, 2, 3, 4, 1}, /* 0 front: minor groove */
      {5, 6, 7, 8, 5}, /* 1 back */
      {1, 4, 8, 5, 1}, /* 2 left */
      {2, 3, 7, 6, 2}, /* 3 right */
      {1, 2, 6, 5, 1}, /* 4 upper */
      {3, 4, 8, 7, 3}  /*  5 lower */
  };
  long default_size[2] = {PS_DFTSIZE, FIG_DFTSIZE}; /* PS & XFIG */
  long is_color, same_faces, updown, mgroove;
  long i, j, k, ioffset8, urxy[3];
  long *depth, *idx, **allobj;
  FILE *fp;

  is_color = opts[2]; /* decomposed for clarity: cf. alc2ps/alc2fig */
  same_faces = opts[7];
  updown = opts[9];
  mgroove = opts[10];

  fp = open_file(imgfile, "w");
  if (opts[0] == 2) { /* Raster3D */
    atom_col = dmatrix(0, NATOMCOL, 1, 3);
    base_col = dmatrix(0, NBASECOL, 1, 3);
    raster3d_header(num_blk8, blkxyz, scale_factor, opts[8], opts[1], fp);
    /* only base_col is used here */
    get_r3dpars(base_col, hb_col, width3, atom_col, label_style);
    for (i = 1; i <= num_blk; i++) {
      ioffset8 = (i - 1) * 8;
      k = is_color ? blkibase[i] : NON_WC_IDX; /* if color image */
      r3d_block_edge(same_faces ? base_col[NON_WC_IDX]
                                : base_col[k],         /* if same-face */
                     ioffset8, blkxyz, width3[2], fp); /* bp1 width */
      for (j = 0; j < 6; j++) {
        if (same_faces)
          pd = (mgroove && !j) ? base_col[NON_WC_IDX] : base_col[k];
        else {
          if ((!updown && !j) || (updown && j == 4))
            pd = base_col[k];
          else
            pd = base_col[NBASECOL];
        }
        r3d_tripln(1, blkxyz[ioffset8 + faces[j][0]],
                   blkxyz[ioffset8 + faces[j][1]],
                   blkxyz[ioffset8 + faces[j][2]], pd, fp); /* 0-1-2 */
        r3d_tripln(1, blkxyz[ioffset8 + faces[j][0]],
                   blkxyz[ioffset8 + faces[j][2]],
                   blkxyz[ioffset8 + faces[j][3]], pd, fp); /* 0-2-3 */
      }
    }
    if (!is_color)
      cpxyz(base_col[NON_WC_IDX], hb_col);
    for (i = 1; i <= nO_lkg; i++) {
      if (!linkage[i][0]) /* bp-center connecting line */
        r3d_dash(oxyz[linkage[i][1]], oxyz[linkage[i][2]], width3[1], hb_col,
                 fp);
      else /* helix or reference frame */
        r3d_rod(3, oxyz[linkage[i][1]], oxyz[linkage[i][2]],
                (linkage[i][0] == -1) ? width3[3] : width3[2], hb_col, fp);
    }

    free_dmatrix(atom_col, 0, NATOMCOL, 1, 3);
    free_dmatrix(base_col, 0, NBASECOL, 1, 3);
  } else {
    nobj += nO_lkg; /* total number of objects */

    /* combine all objects (faces + origin lines) together */
    allobj = lmatrix(1, 3, 1, nobj);
    get_alc_objs(num_blk, blkxyz, nO, oxyz, nO_lkg, linkage, faces, allobj);

    idx = lvector(1, nobj);
    lsort(nobj, allobj[3], idx);

    adjust_xy(num_blk8, blkxyz, nO, oxyz, scale_factor, default_size[opts[0]],
              urxy);

    if (opts[0] == 1) { /* for XFIG */
      depth = lvector(1, nobj);
      get_depth(nobj, allobj[3], depth);
      get_fig_xy(num_blk8, blkxyz, nO, oxyz, urxy, opts[1], fp);
      alc2fig(nobj, idx, depth, allobj, blkxyz, oxyz, blkibase, faces, opts,
              fp);
      free_lvector(depth, 1, nobj);
    } else { /* PS */
      get_ps_xy(imgfile, urxy, opts[1], fp);
      alc2ps(nobj, idx, allobj, blkxyz, oxyz, blkibase, faces, opts, fp);
    }
    free_lmatrix(allobj, 1, 3, 1, nobj);
    free_lvector(idx, 1, nobj);
  }
  close_file(fp);
}

/* "allobj" has three row as follows
 *  (a) for blocks        [1] block index (1 -- num_blk)
 *                        [2] face index (0 -- 5)
 *                        [3] average z-coord x 1000 for each face
 *  (b) for O linkages    [1] 0: origin, -1: helix, -2, -3, -4: frame x, y, z
 *                        [2] 10000 x id1 + id2
 *                        [3] average z-coord x 1000 for the line */
void get_alc_objs(long num_blk, double **blkxyz, long nO, double **oxyz,
                  long nO_lkg, long **linkage, long faces[][5], long **allobj) {
  long i, ioffset8, ip, j, k;

  for (i = 1; i <= num_blk; i++) {
    ioffset8 = (i - 1) * 8; /* ith block atom index offset */
    k = (i - 1) * 6 + 1;    /* ith block face offset */
    for (j = 0; j < 6; j++) {
      ip = k + j;
      allobj[1][ip] = i;
      allobj[2][ip] = j;
      /* only two diagonal atoms are needed for z-coordinates */
      allobj[3][ip] = lround(500.0 * (blkxyz[ioffset8 + faces[j][0]][3] +
                                      blkxyz[ioffset8 + faces[j][2]][3]));
    }
  }

  if (nO) {
    k = num_blk * 6; /* current object offset */
    for (i = 1; i <= nO_lkg; i++) {
      j = k + i;
      allobj[1][j] = linkage[i][0];
      allobj[2][j] = linkage[i][1] * 10000 + linkage[i][2];
      allobj[3][j] =
          lround(500.0 * (oxyz[linkage[i][1]][3] + oxyz[linkage[i][2]][3]));
    }
  }
}

/* read in rendering parameters for XFIG file */
void get_fig_pars(double *dot_sep, long *dlcol, long *dwidth, long *bp1width,
                  long *bp2width, long **bc_idx, double *msat, double *Msat,
                  long *o_sides, long *line_width, long *join_style,
                  long *cap_style, long *mfcol) {
  char BDIR[BUF1K], str[BUF512];
  char *fig_image_par = "fig_image.par";
  long i;
  FILE *fp;

  get_BDIR(BDIR, fig_image_par);
  strcat(BDIR, fig_image_par);
  fprintf(stderr, " ...... reading file: %s ...... \n", fig_image_par);

  fp = open_file(BDIR, "r");
  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", dot_sep) != 1)
    fatal("error in reading dot separation\n");

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", dlcol) != 1)
    fatal("error in reading dot line color\n");

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", dwidth) != 1)
    fatal("error in reading dot line width\n");

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", bp1width) != 1)
    fatal("error in reading bp1 width\n");

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", bp2width) != 1)
    fatal("error in reading bp2 width\n");

  for (i = 0; i <= 6; i++)
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, "%ld", &bc_idx[1][i]) != 1)
      fatal("error in reading color-code\n");
  for (i = 0; i <= 6; i++) /* black & white option */
    bc_idx[0][i] = 0;

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", msat) != 1)
    fatal("error in reading minor groove color saturation\n");
  *msat = fabs(*msat);
  if (*msat > 1.0)
    *msat = 0.9;
  *msat = 1.0 - *msat;
  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", Msat) != 1)
    fatal("error in reading major groove color saturation\n");
  *Msat = fabs(*Msat);
  if (*Msat > 1.0)
    *Msat = 0.1;
  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", o_sides) != 1)
    fatal("error in reading other-side color code\n");

  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", line_width) != 1)
    fatal("error in reading line width\n");
  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", join_style) != 1)
    fatal("error in reading join style\n");
  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", cap_style) != 1)
    fatal("error in reading cap style\n");
  if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", mfcol) != 1)
    fatal("error in reading minor groove filling color\n");

  close_file(fp);
}

/* align the structure so that its global reference frame is in that
 * defined in 'morg' and 'mst' (column-wise for x-, y- and z-axes) */
void frame_xyz(long side_view, double *morg, double **mst, long num,
               double **xyz) {
  double **tmpxyz;
  long i, j;

  tmpxyz = dmatrix(1, num, 1, 3);

  for (i = 1; i <= num; i++)
    for (j = 1; j <= 3; j++)
      tmpxyz[i][j] = dot(xyz[i], mst[j]) + morg[j];

  copy_dmatrix(tmpxyz, num, 3, xyz);

  if (side_view)
    get_side_view(1, num, xyz);

  free_dmatrix(tmpxyz, 1, num, 1, 3);
}

void change_xyz(long side_view, double *morg, double **mst, long num,
                double **xyz) {
  double **tmpxyz;

  tmpxyz = dmatrix(1, num, 1, 3);

  move_position(xyz, num, 3, morg);
  multi_matrix(xyz, num, 3, mst, 3, 3, tmpxyz);

  copy_dmatrix(tmpxyz, num, 3, xyz);

  if (side_view)
    get_side_view(1, num, xyz);

  free_dmatrix(tmpxyz, 1, num, 1, 3);
}

/* adjust orientation: xyz*rotz(-90)*rotx(90) ==> [-y z -x] */
void get_side_view(long ib, long ie, double **xyz) {
  long i;
  double temp;

  for (i = ib; i <= ie; i++) {
    temp = xyz[i][1];
    xyz[i][1] = -xyz[i][2];
    xyz[i][2] = xyz[i][3];
    xyz[i][3] = -temp;
  }
}

/* get C1* and RN9/YN1 atom pair index for finding helix */
void get_CNidx(long ds, long num_bp, long **chi, long **idx, long *nvec,
               long *C1b, long *C1e) {
  long i, ia, ib, ioffset, j, joffset, k;

  *nvec = 0;
  for (i = 1; i <= ds; i++)
    for (j = 1; j <= num_bp - 1; j++) {
      ioffset = (j - 1) * 4;
      joffset = j * 4;
      for (k = 2; k <= 3; k++) {
        ia = chi[i][ioffset + k];
        ib = chi[i][joffset + k];
        if (ia && ib) {
          idx[++*nvec][1] = ia;
          idx[*nvec][2] = ib;
        }
      }
    }

  ioffset = (num_bp - 1) * 4;
  for (i = 1; i <= ds; i++) {
    C1b[i] = chi[i][2];
    C1e[i] = chi[i][ioffset + 2];
  }
}

/* attach reference frame axes */
void add_3axes(long *num, char **AtomName, long *ibase, double **xyz,
               long *nbond, long **linkage, long side_view, double axis_len) {
  static char *ref_symbol[4] = {"OO", "OX", "OY", "OZ"};
  static double ref_xyz[4][3] = {
      {0.0, 0.0, 0.0}, /* origin */
      {1.0, 0.0, 0.0}, /* x-axis */
      {0.0, 1.0, 0.0}, /* y-axis */
      {0.0, 0.0, 1.0}  /* z-axis */
  };
  long i, j, k;

  for (i = 1; i <= 4; i++) {
    k = *num + i;
    strcpy(AtomName[k], ref_symbol[i - 1]);
    for (j = 1; j <= 3; j++)
      xyz[k][j] = axis_len * ref_xyz[i - 1][j - 1];
    ibase[k] = NON_WC_IDX; /* non-common base atoms */
  }
  if (side_view)
    get_side_view(*num + 1, *num + 4, xyz);
  for (i = 1; i <= 3; i++) {
    k = *nbond + i;
    linkage[k][1] = *num + 1;
    linkage[k][2] = *num + i + 1;
  }
  *num += 4;
  *nbond += 3;
}

char *get_sequence(char *Wbase, long *num_bp) {
  char *bseq;
  char str[BUF512];
  long ich = 2, k, nbp;

  fprintf(stderr, "\nInput your base sequence with only %s:\n", Wbase);
  fprintf(stderr, "1. From a data file (complete sequence)\n");
  fprintf(stderr, "2. From keyboard (enter only the repeating sequence)\n");
  fprintf(stderr, "Your choice (1 or 2, Dft: 2): ");
  fflush(stderr);
  if (fgets(str, sizeof str, stdin) != NULL) {
    k = sscanf(str, "%ld", &ich);
    if (!k || k == EOF)
      ich = 2;
  } else
    fatal("error in reading your choice\n");
  fprintf(stderr, "\n");

  if (ich == 1) {
    fprintf(stderr, "Name of your base sequence file: ");
    fflush(stderr);
    if (fgets(str, sizeof str, stdin) == NULL)
      fatal("error in reading your sequence file name\n");
    trim(str);
    bseq = read_sequence(str, Wbase, &nbp);
  } else
    bseq = read_repeat("A", 0, Wbase, &nbp);

  *num_bp = nbp;

  return bseq;
}

char **single2double(long nbp, char *bseq, char *Wbase, char *Cbase) {
  char **bp_seq;
  long i, j;

  bp_seq = cmatrix(1, nbp, 1, 2);
  for (i = 1; i <= nbp; i++) {
    bp_seq[i][1] = bseq[i - 1];
    j = strchr(Wbase, bseq[i - 1]) - Wbase;
    bp_seq[i][2] = *(Cbase + j);
  }
  free_cvector(bseq, 0, nbp);

  return bp_seq;
}

long is_valid_base(char c, char *valid_bases) {
  if (isspace((int)c)) /* skip white space */
    return false;

  if (strchr(valid_bases, c) == NULL) {
    fprintf(stderr, "skip %c: acceptable bases [%s]\n", c, valid_bases);
    return false;
  }

  return true;
}

long repeat_num(void) {
  char *p0;
  long num = 10, k;

  fprintf(stderr, "Number of repeats (Dft: 10): ");
  fflush(stderr);
  if ((p0 = my_getline(stdin)) != NULL) {
    k = sscanf(p0, "%ld", &num);
    if (!k || k == EOF || !num)
      num = 10;
  } else
    fatal("error in reading number of repeats\n");
  free(p0);

  if (num < 0) {
    fprintf(stderr, "number of repeat %ld < 0. reset to positive\n", num);
    num = -num;
  }

  return num;
}

/* read in sequence information from file 'seqfile', return number of
 * valid bases [must be in UPPER case], and array of the bases */
char *read_sequence(char *seqfile, char *valid_bases, long *nbp) {
  char *p0, *line, *pb, *bseq;
  long nb = 0, maxline = BUF512;
  FILE *fp;

  pb = (char *)malloc(maxline * sizeof(char)); /* initial size */
  if (pb == NULL)
    fatal("malloc failure in reading sequence\n");

  fp = open_file(seqfile, "r");
  while ((p0 = my_getline(fp)) != NULL) {
    line = trim(p0); /* keep the original value of p0 */
    if (line[0] == '#' || line[0] == '\0') {
      free(p0);
      continue; /* skip empty and commented lines */
    }
    upperstr(line); /* changed to upper case */
    while (*line) {
      if (is_valid_base(*line, valid_bases)) {
        if (nb >= maxline - 1)
          pb = enlarge_cline(&maxline, pb);
        pb[nb++] = *line;
      }
      line++;
    }
    free(p0);
  }
  pb[nb] = '\0';

  close_file(fp);
  bseq = my_strdup(pb);
  free(pb);

  if (!nb)
    fatal("sequence file <%s> contains no valid bases\n", seqfile);

  *nbp = nb;

  return bseq;
}

/* read in repeated sequence, return number of valid bases [must be in
 * UPPER case], and array of the bases */
char *read_repeat(char *crepeat, long fixed, char *valid_bases, long *nbp) {
  char *p0, *line, *pb, *bseq;
  long i, k, num, num_total, nb = 0, maxline = BUF512;

  pb = (char *)malloc(maxline * sizeof(char));
  if (pb == NULL)
    fatal("malloc failure in reading sequence\n");

  if (fixed) {
    fprintf(stderr, "Repeating unit: %s\n", crepeat);
    fflush(stderr);
    strcpy(pb, crepeat);
    nb = strlen(crepeat);
  } else {
    fprintf(stderr, "Repeating unit (Dft: %s): ", crepeat);
    fflush(stderr);
    if ((p0 = my_getline(stdin)) == NULL)
      fatal("error in reading your repeating unit\n");
    line = trim(p0);
    upperstr(line);

    while (*line) {
      if (is_valid_base(*line, valid_bases)) {
        if (nb >= maxline - 1)
          pb = enlarge_cline(&maxline, pb);
        pb[nb++] = *line;
      }
      line++;
    }
    free(p0);

    if (!nb) { /* using default */
      strcpy(pb, crepeat);
      nb = strlen(crepeat);
    } else
      pb[nb] = '\0';

    fprintf(stderr, "Repeating unit: %s\n", pb);
  }

  num = repeat_num();

  num_total = nb * num;
  bseq = cvector(0, num_total);

  for (i = 1; i <= num; i++) {
    k = (i - 1) * nb;
    strcpy(bseq + k, pb);
  }

  free(pb);

  *nbp = num_total;

  return bseq;
}

/* combine strand II with strand I for a parallel duplex */
void combine_pstnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                    char **tAtomName, char **tResName, char *tChainID,
                    long *tResSeq, double **txyz, char **tAtomName2,
                    double **txyz2) {
  long i, ik, j;

  for (i = 1; i <= num_bp; i++) {
    for (j = s2idx[i][1] + 1; j <= s2idx[i][1] + s2idx[i][2]; j++) {
      ik = *tnum + j - s2idx[i][1];
      strcpy(tAtomName[ik], tAtomName2[j]);
      sprintf(tResName[ik], "  %c", bp_seq[i][2]);
      tChainID[ik] = Gvars.REBUILD_CHAIN_IDS[1];
      tResSeq[ik] = num_bp + i;
      cpxyz(txyz2[j], txyz[ik]);
    }
    *tnum += s2idx[i][2];
  }
}

/* reverse strand II to be in 5'-->3' direction & combined with I */
void reverse_stnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                   char **tAtomName, char **tResName, char *tChainID,
                   long *tResSeq, double **txyz, char **tAtomName2,
                   double **txyz2, long basep) {
  long i, ik, j;

  for (i = num_bp; i >= 1; i--) {
    for (j = s2idx[i][1] + 1; j <= s2idx[i][1] + s2idx[i][2]; j++) {
      ik = *tnum + j - s2idx[i][1];
      strcpy(tAtomName[ik], tAtomName2[j]);
      sprintf(tResName[ik], "  %c", bp_seq[i][2]);
      tChainID[ik] = Gvars.REBUILD_CHAIN_IDS[1];
      tResSeq[ik] = 2 * num_bp - i + 1;
      if (basep && !strcmp(tAtomName[ik], " P  "))
        ++tResSeq[ik]; /* re-position P on strand II */
      cpxyz(txyz2[j], txyz[ik]);
    }
    *tnum += s2idx[i][2];
  }
}

void pair_checking(long ip, long ds, long num_residue, char *pdbfile,
                   long *num_bp, long **pair_num) {
  long i, j;

  if (!ip) { /* ideal base-pairing */
    if (ds == 1) {
      if (*num_bp > num_residue) {
        *num_bp = num_residue;
        fprintf(stderr, "processing all %ld residues\n\n", *num_bp);
      }
    } else if (num_residue % 2) { /* ds = 2 */
      fprintf(stderr, "%s has odd number (%ld) of residues\n", pdbfile,
              num_residue);
      fatal("Please specify base-pairing residue numbers\n");
    } else {
      i = num_residue / 2;
      if (*num_bp < i)
        fatal("Please specify base-pairing residue numbers\n");
      if (*num_bp > i) {
        *num_bp = i;
        fprintf(stderr, "processing all %ld base-pairs\n\n", i);
      }
    }

    for (i = 1; i <= *num_bp; i++) {
      pair_num[1][i] = i;
      if (ds == 2)
        pair_num[2][i] = 2 * *num_bp - i + 1;
    }
  } else {
    if (ds * *num_bp > num_residue)
      fprintf(stderr, "some residue has more than one pair\n");
    for (i = 1; i <= ds; i++)
      for (j = 1; j <= *num_bp; j++)
        if (pair_num[i][j] > num_residue) {
          fprintf(stderr, "residue index %ld too big (> %ld)\n", pair_num[i][j],
                  num_residue);
          fatal("please check your input file\n");
        }
  }
}

void double_print_msg(char *msg, FILE *fp) {
  fprintf(stderr, "%s\n", msg);
  fprintf(fp, "%s\n", msg);
}

/* check for proper O3'-P linkage, strand I 5'--> 3' direction, parallel duplex
 */
void drct_checking(long ds, long num_bp, long **pair_num, long **seidx,
                   char **AtomName, double **xyz, long *parallel, long *bbexist,
                   long **o3p_brk, FILE *fp) {
  char msg[BUF512];
  long i, ib, ie, ioffset, j, jm1, jp1, rnum, n = 0;
  long direction[7], **O3, **P;

  O3 = lmatrix(1, 2, 1, num_bp);
  P = lmatrix(1, 2, 1, num_bp);

  for (i = 1; i <= ds; i++)
    for (j = 1; j <= num_bp; j++) {
      rnum = pair_num[i][j];
      ib = seidx[rnum][1];
      ie = seidx[rnum][2];
      P[i][j] = find_1st_atom(" P  ", AtomName, ib, ie, "");
      O3[i][j] = find_1st_atom(" O3'", AtomName, ib, ie, "");
      if (P[i][j] && O3[i][j]) {
        ++*bbexist;
        if (within_limits(xyz[O3[i][j]], xyz[P[i][j]], 0.8, O3P_UPPER))
          n++;
      }
    }

  if (n) {
    print_sep(fp, '*', 76);
    sprintf(msg,
            "WARNING: %ld out of %ld bases have O3'[i] wrongly connected"
            " to P[i]",
            n, ds * num_bp);
    double_print_msg(msg, fp);
  }

  init_lvector(direction, 1, 6, 0);
  for (i = 1; i <= ds; i++) {
    ioffset = (i - 1) * 3;
    for (j = 1; j < num_bp; j++) {
      jp1 = j + 1;
      if (O3[i][j] && P[i][jp1] && /* j to j + 1 */
          within_limits(xyz[P[i][jp1]], xyz[O3[i][j]], 0.8, O3P_UPPER)) {
        ++direction[ioffset + 1];
        continue;
      }
      if (O3[i][jp1] && P[i][j] && /* j + 1 to j */
          within_limits(xyz[P[i][j]], xyz[O3[i][jp1]], 0.8, O3P_UPPER)) {
        ++direction[ioffset + 2];
        continue;
      }
      ++direction[ioffset + 3];
      if (*bbexist)
        o3p_brk[i][j] = 9;
    }
  }

  if ((direction[1] && direction[2]) || (direction[4] && direction[5])) {
    sprintf(msg, "This structure contains intra-chain direction reverse");
    double_print_msg(msg, fp);
  } else {
    if (*bbexist && (direction[3] || direction[6])) {
      sprintf(msg, "This structure has broken O3'[i] to P[i+1] linkages");
      double_print_msg(msg, fp);
    }
    if (direction[2] && !direction[1]) {
      sprintf(msg, "WARNING: Strand I in 3'-->5' direction!!!\n");
      double_print_msg(msg, fp);
    }
    if (direction[4] && !direction[5]) {
      fprintf(stderr, "This is a parallel duplex\n");
      *parallel = 1;
    }
  }
  print_sep(fp, '*', 76);

  if (ds == 1) {
    long ka, kb, ksum;

    for (j = 1; j <= num_bp; j++) {
      ksum = 0;

      jm1 = j - 1;
      if (jm1 < 1)
        ksum++;
      else {
        ka = O3[1][jm1];
        kb = P[1][j];
        if (ka && kb && within_limits(xyz[ka], xyz[kb], 0.8, O3P_UPPER))
          continue;
        ksum++;
      }

      jp1 = j + 1;
      if (jp1 > num_bp)
        ksum++;
      else {
        ka = O3[1][j];
        kb = P[1][jp1];
        if (ka && kb && within_limits(xyz[ka], xyz[kb], 0.8, O3P_UPPER))
          continue;
        ksum++;
      }

      if (*bbexist && ksum == 2)
        o3p_brk[1][j] = 1;
    }
  }

  free_lmatrix(O3, 1, ds, 1, num_bp);
  free_lmatrix(P, 1, ds, 1, num_bp);
}

void residue_idstr(char chain_id, long res_seq, char *rname, char *idmsg) {
  long i;

  sprintf(idmsg, "%c:%ld:%s", chain_id, res_seq, rname);
  for (i = 0; i < (long)strlen(idmsg); i++)
    if (idmsg[i] == ' ')
      idmsg[i] = '_';
}

/* get base name: all information to uniquely identify a base residue */
void base_str(char chain_id, long res_seq, char *misc, char *rname, char bcode,
              long stnd, char *idmsg) {
  char b1, snum[10], rname_cp[10], modelNum[10], iCode = misc[2];
  long i, nlen = 4;

  b1 = (bcode == '\0') ? ' ' : bcode;

  strncpy(modelNum, misc + 30, nlen);
  modelNum[nlen] = '\0';
  for (i = 0; i < nlen; i++)
    if (modelNum[i] == '0')
      modelNum[i] = '.';
  sprintf(snum, "%4ld", res_seq);
  for (i = 0; i < 4; i++)
    if (snum[i] == ' ')
      snum[i] = '.';
  if (chain_id == ' ')
    chain_id = '-';
  if (iCode == ' ')
    iCode = '_';
  strcpy(rname_cp, rname);
  for (i = 0; i < 3; i++)
    if (rname_cp[i] == ' ')
      rname_cp[i] = '.';
  if (stnd == 1) /* strand I */
    sprintf(idmsg, "%4s>%c:%4s%c:[%s]%c", modelNum, chain_id, snum, iCode,
            rname_cp, b1);
  else /* strand II */
    sprintf(idmsg, "%c[%s]:%4s%c:%c<%4s", b1, rname_cp, snum, iCode, chain_id,
            modelNum);
}

/* write out a bond list for checking and easy modification */
void write_lkglist(long nbond, long **linkage, char **AtomName, char **ResName,
                   char *ChainID, long *ResSeq, char **Miscs) {
  char b1[BUF512], b2[BUF512];
  long i, ia, ib;
  FILE *fp;

  fp = open_file(LKG_FILE, "w");
  for (i = 1; i <= nbond; i++) {
    ia = linkage[i][1];
    ib = linkage[i][2];
    base_str(ChainID[ia], ResSeq[ia], Miscs[ia], ResName[ia], '\0', 1, b1);
    base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], '\0', 1, b2);
    fprintf(fp, "%5ld %5ld  # %5ld %s {%s} with %s {%s}\n", ia, ib, i,
            AtomName[ia], b1, AtomName[ib], b2);
  }
  close_file(fp);
}

/* get all possible H-bonds with upper/lower layer information */
void hbond_info(long **pair_num, char *bseq, long **seidx, long *idx,
                char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, long *RY, long *num_hbond,
                long **hb_linkage) {
  long i, j, k, ki, kj, num_bp = 2;
  long **num_list;
  FILE *fph;

  num_list = lmatrix(1, BUF512, 0, 3);
  fph = open_file(HB_FILE, "w");
  for (k = 1; k <= num_bp; k++) {
    for (i = 1; i < pair_num[k][0]; i++) {
      ki = pair_num[k][i];
      if (RY[ki] < 0)
        continue;
      for (j = i + 1; j <= pair_num[k][0]; j++) { /* pairwise */
        kj = pair_num[k][j];
        if (RY[kj] < 0)
          continue;
        hbond_list(ki, kj, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, bseq,
                   seidx, idx, hb_linkage, &Gvars.misc_pars, num_list,
                   num_hbond, k, fph);
      }
    }
  }
  close_file(fph);
  free_lmatrix(num_list, 1, BUF512, 0, 3);
}

/* get H-bonds betwen base-pairs for a general PDB data file */
void hbond_pdb(long num, long num_residue, char *bseq, long **seidx, long *idx,
               long *RY, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, double **xyz, long *num_hbond,
               long **hb_linkage, long pwise) {
  double rtn_val[RTNNUM];
  double **orien, **org, **NC1xyz, **o3_p;
  long bpid, i, j, **num_list, **ring_atom;
  FILE *fph;

  orien = dmatrix(1, num_residue, 1, 9);
  org = dmatrix(1, num_residue, 1, 3);
  NC1xyz = dmatrix(1, num_residue, 1, 7); /* RN9/YN1 & C1' atomic coordinates */
  o3_p = dmatrix(1, num_residue, 1, 8);   /* O3'/P atomic coordinates */

  base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq,
            Miscs, xyz, orien, org, NC1xyz, o3_p);

  /* 1-9 ring atom index, 10 # of ring atoms, 11-19 first level */
  ring_atom = lmatrix(1, num_residue, 1, 19);
  ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);

  num_list = lmatrix(1, BUF512, 0, 3);

  fph = open_file(HB_FILE, "w");
  for (i = 1; i < num_residue; i++) {
    if (!pwise && RY[i] < 0)
      continue;
    for (j = i + 1; j <= num_residue; j++) {
      if (!pwise) {    /* check for base pairing */
        if (RY[j] < 0) /* only considering hydration for the 1st layer */
          continue;
        check_pair(i, j, bseq, seidx, xyz, NC1xyz, orien, org, idx, AtomName,
                   &Gvars.misc_pars, rtn_val, &bpid, ring_atom, 0);
        if (!bpid)
          continue; /* not a base-pair */
      }
      hbond_list(i, j, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, bseq,
                 seidx, idx, hb_linkage, &Gvars.misc_pars, num_list, num_hbond,
                 1, fph);
    }
  }

  free_dmatrix(orien, 1, num_residue, 1, 9);
  free_dmatrix(org, 1, num_residue, 1, 3);
  free_dmatrix(NC1xyz, 1, num_residue, 1, 7);
  free_dmatrix(o3_p, 1, num_residue, 1, 8);
  free_lmatrix(ring_atom, 1, num_residue, 1, 19);
  free_lmatrix(num_list, 1, BUF512, 0, 3);
  close_file(fph);
}

/* write out h-bond listing information: <num_hbond> already a pointer */
void hbond_list(long i, long j, char **AtomName, char **ResName, char *ChainID,
                long *ResSeq, double **xyz, char **Miscs, char *bseq,
                long **seidx, long *idx, long **hb_linkage, miscPars *misc_pars,
                long **num_list, long *num_hbond, long ilayer, FILE *fph) {
  char b1[BUF512], b2[BUF512];
  long ir, jr, k, num_hb;

  ir = seidx[i][1];
  jr = seidx[j][1];
  base_str(ChainID[ir], ResSeq[ir], Miscs[ir], ResName[ir], bseq[i], 1, b1);
  base_str(ChainID[jr], ResSeq[jr], Miscs[jr], ResName[jr], bseq[j], 1, b2);
  hb_numlist(i, j, bseq[i], bseq[j], seidx, idx, AtomName, xyz, misc_pars,
             &num_hb, num_list);
  for (k = 1; k <= num_hb; k++) {
    hb_linkage[++*num_hbond][1] = ilayer;
    hb_linkage[*num_hbond][2] = num_list[k][1];
    hb_linkage[*num_hbond][3] = num_list[k][2];
    fprintf(fph, "%5ld %5ld  # %5ld %4.2f%c", num_list[k][1], num_list[k][2],
            *num_hbond, num_list[k][3] / MFACTOR, num_list[k][0] ? '*' : ' ');
    fprintf(fph, "%s {%s} with %s {%s}\n", AtomName[num_list[k][1]], b1,
            AtomName[num_list[k][2]], b2);
  }
}

/* determine atom number list for each base-pair from string hb_info */
void hb_numlist(long i, long j, char basei, char basej, long **seidx, long *idx,
                char **AtomName, double **xyz, miscPars *misc_pars,
                long *num_hb, long **num_list) {
  char atom1[10], atom2[10], hb_info[BUF512], *pchar;
  long k, num;

  get_hbond_ij(i, j, basei, basej, misc_pars, seidx, idx, AtomName, xyz,
               hb_info);

  if (sscanf(hb_info, "[%ld]", &num) != 1)
    num = 0;
  if (num) {
    pchar = strchr(hb_info, ' ');
    for (k = 1; k <= num; k++) {
      strncpy(atom1, pchar + 1, 4);
      atom1[4] = '\0';
      strncpy(atom2, pchar + 6, 4);
      atom2[4] = '\0';
      num_list[k][1] =
          find_1st_atom(atom1, AtomName, seidx[i][1], seidx[i][2], "");
      num_list[k][2] =
          find_1st_atom(atom2, AtomName, seidx[j][1], seidx[j][2], "");

      num_list[k][0] = 0;
      if (*(pchar + 5) == '*')
        num_list[k][0] = 1;          /* not a canonical D-A atom pair */
      strncpy(atom1, pchar + 11, 4); /* HB distance */
      atom1[4] = '\0';
      num_list[k][3] = lround(cvt2double(atom1) * MFACTOR);
      pchar += 15; /* point to next H-bond */
    }
  }
  *num_hb = num;
}

/* get detailed H-bonding information for each pair */
void hb_information(long num_bp, long **pair_num, char **bp_seq, long **seidx,
                    long *idx, char **AtomName, double **xyz, long *WC_info,
                    FILE *fp) {
  char hb_info[BUF512];
  long i, j, k;

  print_sep(fp, '*', 76);
  fprintf(fp, "Detailed H-bond information: atom-name pair and length [%s]\n",
          Gvars.misc_pars.hb_atoms);

  for (k = 1; k <= num_bp; k++) {
    i = pair_num[1][k];
    j = pair_num[2][k];
    get_hbond_ij(i, j, bp_seq[1][k], bp_seq[2][k], &Gvars.misc_pars, seidx, idx,
                 AtomName, xyz, hb_info);
    fprintf(fp, "%4ld %c-%c%c%c-%c  %s\n", k, bp_seq[1][k],
            (WC_info[k] == 2) ? '-' : '*', WC_info[k] ? '-' : '*', bp_seq[0][k],
            bp_seq[2][k], hb_info);
  }
}

/* check if two atom pair could form a H-bond */
long good_hbatoms(miscPars *misc_pars, char *atom1, char *atom2, long idx1,
                  long idx2) {
  static char *PO[] = {" O1P", " O2P", " O3'", " O4'", " O5'", " N7 "};
  static long numPO = sizeof PO / sizeof PO[0] - 1;
  long natom = misc_pars->hb_idx[0];

  if (num_strmatch(atom1, PO, 0, numPO) && num_strmatch(atom2, PO, 0, numPO))
    return false; /* no H-bond within PO4 group, O4' & N7 */

  if ((idx1 == 2 || idx1 == 4 || idx2 == 2 ||
       idx2 == 4) && /* at least one O/N atom */
      (lval_in_set(idx1, 1, natom, misc_pars->hb_idx) &&
       lval_in_set(idx2, 1, natom, misc_pars->hb_idx)))
    return true;
  else
    return false;
}

/* read in bond linkage information from "lkgfile": r3d_atom & stack2img
 * generate a file called "bonds_lkg.dat", which can be modified to fit any
 * purpose */
void read_lkginfo(char *lkgfile, long num, long *nbond, long **linkage) {
  char str[BUF512];
  long ia, ib;
  FILE *fp;

  fp = open_file(lkgfile, "r");
  while (fgets(str, sizeof str, fp) != NULL)
    if (sscanf(str, "%ld %ld", &ia, &ib) == 2) {
      if (!lval_in_range(ia, 1, num) || !lval_in_range(ib, 1, num)) {
        fprintf(stderr, "%s", str);
        fprintf(stderr, "         has atom serial number out of range\n");
        continue;
      }
      linkage[++*nbond][1] = ia;
      linkage[*nbond][2] = ib;
    }
  close_file(fp);
}

/* read in H-bonding information from "hbfile": r3d_atom & stack2img generate
 * a file called "hbonds_info.dat", which can be modified to fit any purpose */
void read_hbinfo(char *hbfile, long num, long *num_hbond, long **hb_linkage) {
  char str[BUF512];
  long ia, ib;
  FILE *fp;

  fp = open_file(hbfile, "r");
  while (fgets(str, sizeof str, fp) != NULL)
    if (sscanf(str, "%ld %ld", &ia, &ib) == 2) {
      if (!lval_in_range(ia, 1, num) || !lval_in_range(ib, 1, num)) {
        fprintf(stderr, "[%s] has atom serial number out of range\n", str);
        continue;
      }
      hb_linkage[++*num_hbond][1] = 1;
      hb_linkage[*num_hbond][2] = ia;
      hb_linkage[*num_hbond][3] = ib;
    }
  close_file(fp);
}

void update_hb_idx(long idx, double *dtmp, long *ddidx, double *hb_dist,
                   long cur_idx) {
  dtmp[idx] = hb_dist[cur_idx];
  ddidx[idx] = cur_idx;
}

/* get the H-bonding atom-pairs */
void hb_atompair(long num_hbonds, char **hb_atom1, char **hb_atom2,
                 double *hb_dist, long *lkg_type, miscPars *misc_pars) {
  double dtmp[3];
  long k, m = 0, n, num_iter = 1;
  long ddidx[3], *matched_idx, **idx2;

  if (!num_hbonds)
    return;

  matched_idx = lvector(1, num_hbonds);

  while (1) {
    if (matched_idx[num_iter]) {
      num_iter++;
      continue;
    }

    for (k = 1; k <= 2; k++)
      update_hb_idx(k, dtmp, ddidx, hb_dist, num_iter);

    for (n = 1; n <= num_hbonds; n++) {
      if (n == num_iter || matched_idx[n])
        continue;
      if (!strcmp(hb_atom1[n], hb_atom1[num_iter]) && hb_dist[n] < dtmp[1])
        update_hb_idx(1, dtmp, ddidx, hb_dist, n);
      if (!strcmp(hb_atom2[n], hb_atom2[num_iter]) && hb_dist[n] < dtmp[2])
        update_hb_idx(2, dtmp, ddidx, hb_dist, n);
    }

    if (ddidx[1] == ddidx[2]) { /* best mutual match */
      k = ddidx[1];
      hb_dist[k] = -hb_dist[k]; /* make  */
      num_iter = 1;             /* reset the iterator */

      for (n = 1; n <= num_hbonds; n++) {
        if (matched_idx[n])
          continue;
        if (!strcmp(hb_atom1[n], hb_atom1[k]) || /* not && */
            !strcmp(hb_atom2[n], hb_atom2[k])) {
          matched_idx[n] = 1;
          m++;
        }
      }

      if (m >= num_hbonds)
        break;

    } else
      num_iter++;
  }

  /* === further processing by adding more possible H-bonds === */
  idx2 = lmatrix(1, num_hbonds, 1, 2);

  for (k = 1; k <= num_hbonds; k++) {
    if (hb_dist[k] > 0.0)
      continue;
    idx2[k][1] = 9;
    idx2[k][2] = 9;
    for (m = 1; m <= num_hbonds; m++) {
      if (m == k || hb_dist[m] < 0.0)
        continue;
      if (!strcmp(hb_atom1[m], hb_atom1[k]))
        idx2[m][1] = 1;
      if (!strcmp(hb_atom2[m], hb_atom2[k]))
        idx2[m][2] = 1;
    }
  }

  /* Note that at this point, there are only 3 possibilities:
     9 + 9 = 18 means already located H-bonds;
     1 + 1 = 2 means the two atoms in separate H-bonds;
     1 + 0 = 0 + 1 = 1 means only one of atoms in H-bonds */
  for (k = 1; k <= num_hbonds; k++) {
    m = idx2[k][1] + idx2[k][2];
    lkg_type[k] = m; /* for later on referencing */
    if (m != 18 &&
        dval_in_range(hb_dist[k], misc_pars->hb_lower, misc_pars->hb_dist2))
      hb_dist[k] = -hb_dist[k];
  }

  /* The following code for adding additional H-bonds is too
   * complicated, ad hoc, but seems adding not much useful information
   * as from the above simple approach -- feb28, 2008 */

/* ============================================================================
 */
#if 0
    for (k = 1; k <= num_hbonds; k++) {
        m = idx2[k][1] + idx2[k][2];
        lkg_type[k] = m;  /* for later on referencing */
        if (m == 1 && dval_in_range(hb_dist[k], misc_pars->hb_lower, misc_pars->hb_dist2))
            hb_dist[k] = -hb_dist[k];
    }

    /* post processing to account for cases such as (I34 T--A J259) 1zla -- nov-02-2006 */
    for (m = 1; m <= num_hbonds; m++) {
        matched_idx[m] = 0;
        k = idx2[m][1] + idx2[m][2];
        if (k != 18)
            continue;
        for (n = 1; n <= num_hbonds; n++) {
            if (hb_dist[n] < 0.0)  /* covering n == m, 9 + 9 == 18, 0 + 1 etc */
                continue;
            if (hb_dist[n] > fabs(hb_dist[m]) ||
                !dval_in_range(hb_dist[n], misc_pars->hb_lower, misc_pars->hb_dist2))
                continue;
            if (!strcmp(hb_atom1[m], hb_atom1[n]) || !strcmp(hb_atom2[m], hb_atom2[n])) {
                hb_dist[n] = -hb_dist[n];
                matched_idx[n] = 1;
            }
        }
    }

    for (k = 1; k <= num_hbonds; k++)
        if (matched_idx[k] && !dval_in_range(fabs(hb_dist[k]), misc_pars->hb_lower,
                                             misc_pars->hb_dist2))
            hb_dist[k] = fabs(hb_dist[k]);
#endif
  /* ============================================================================
   */

  free_lmatrix(idx2, 1, num_hbonds, 1, 2);
  free_lvector(matched_idx, 1, num_hbonds);

  /* === End of further processing by adding more possible H-bonds === */
}

long validate_hbonds(long num_hbonds, double *hb_dist, long *lkg_type,
                     char *hb_type, char basei, char basej, char **hb_atom1,
                     char **hb_atom2) {
  long k, m = 0;

  for (k = 1; k <= num_hbonds; k++) {
    hb_type[k] = ' ';     /* initialized to ' ' */
    if (hb_dist[k] > 0.0) /* not taken as H-bond */
      continue;
    hb_type[k] = donor_acceptor(basei, basej, hb_atom1[k], hb_atom2[k]);
    hb_dist[k] = fabs(hb_dist[k]); /* using hb_type as an indicator */
    if (hb_type[k] == '-' && dval_in_range(hb_dist[k], 2.5, 3.5))
      m++;
  }

  if (m) { /* with good D-A type H-bond, further checking & elimination */
    for (k = 1; k <= num_hbonds; k++) {
      if (hb_type[k] == ' ')
        continue;
      if ((hb_dist[k] > 3.6) || /* must be within reasonable limit 3.6 */
          (hb_type[k] == '*' &&
           lkg_type[k] != 18 && /* D-D or A-A, not in best link */
           !dval_in_range(hb_dist[k], 2.6,
                          3.2))) /* e.g. A1498[U-**+-A]A1499 of 1vs7 */
        hb_type[k] = ' ';
    }
  }

  /* now count and output 'good' H-bonding information */
  m = 0;
  for (k = 1; k <= num_hbonds; k++) {
    if (hb_type[k] == ' ')
      continue;
    m++;
  }

  return m;
}

/* get H-bond length information between residue i and j */
void get_hbond_ij(long i, long j, char basei, char basej, miscPars *misc_pars,
                  long **seidx, long *idx, char **AtomName, double **xyz,
                  char *hb_info) {
  char *hb_type, **hb_atom1, **hb_atom2, aname1[5], aname2[5], stmp[20];
  double *hb_dist;
  long k, m, n, num_hbonds = 0, *lkg_type;

  hb_atom1 = cmatrix(1, BUF512, 0, 4);
  hb_atom2 = cmatrix(1, BUF512, 0, 4);
  hb_dist = dvector(1, BUF512);

  for (m = seidx[i][1]; m <= seidx[i][2]; m++) {
    for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
      if (good_hbatoms(misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]) &&
          within_limits(xyz[n], xyz[m], misc_pars->hb_lower,
                        misc_pars->hb_dist1)) {
        if (++num_hbonds > BUF512)
          fatal("Too many possible H-bonds between two bases\n");
        strcpy(hb_atom1[num_hbonds], AtomName[m]);
        strcpy(hb_atom2[num_hbonds], AtomName[n]);
        hb_dist[num_hbonds] = p1p2_dist(xyz[n], xyz[m]);
      }
    }
  }

  /* Yurong reported a case in rr0027 (1i97):
     >A:..39_:[..G]G-**+-A[..A]:.530_:A< where the bp criteria defines it as a
     pair, because G:N7-A:N7 distance is 3.89A which is within 4.0A cut-off.
     However, in this function when H-bonding info is checked, there is no
     good_hbatoms() pairs. Thus num_hbonds is ZERO and hb_info[] would be
     undefined, leaving it to whatever content the system assign it with --
     mostly garbage. Here we explicitly initialize it to an empty string to
     avoid this problem. */
  if (!num_hbonds) /* no H-bond found */
    sprintf(hb_info, "[%ld]", num_hbonds);

  else {
    lkg_type = lvector(1, num_hbonds);
    hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, misc_pars);

    /*  to eliminate erroneous H-bonds, based purely on geometrical criteria */
    hb_type = cvector(1, num_hbonds);
    m = validate_hbonds(num_hbonds, hb_dist, lkg_type, hb_type, basei, basej,
                        hb_atom1, hb_atom2);

    sprintf(hb_info, "[%ld]", m);
    for (k = 1; k <= num_hbonds; k++) {
      if (hb_type[k] == ' ')
        continue;
      strcpy(aname1, hb_atom1[k]);
      cvt_pdbv3_name(aname1);
      strcpy(aname2, hb_atom2[k]);
      cvt_pdbv3_name(aname2);
      sprintf(stmp, " %s%c%s %4.2f", aname1, hb_type[k], aname2, hb_dist[k]);
      strcat(hb_info, stmp);
    }

    free_cvector(hb_type, 1, num_hbonds);
    free_lvector(lkg_type, 1, num_hbonds);
  }

  free_cmatrix(hb_atom1, 1, BUF512, 0, 4);
  free_cmatrix(hb_atom2, 1, BUF512, 0, 4);
  free_dvector(hb_dist, 1, BUF512);
}

/* based on canonical form of ACGITU, decide if a HB-pair is donor-acceptor */
char donor_acceptor(char basei, char basej, char *hb_atom1, char *hb_atom2) {
  char da[3], hbatom_type = '*', *pchar;
  char ia = '\0', ja = '\0';
  /* _A: acceptor; _D: donor; _X: both acceptor & donor; _?: not sure */
  static char *cmn_base = CB_LIST;
  static char *da_type[7] = {"AD", "AX", "XD", "XX", "DA", "DX", "XA"};
  static char *bb_da[6] = {" O1P_A", " O2P_A", " O5'_A",
                           " O4'_A", " O3'_A", " O2'_X"};
  static char *base_da[6][6] = {
      {" N9 _?", " N7 _A", " N6 _D", " N1 _A", " N3 _A"},           /* A */
      {" N1 _?", " O2 _A", " N3 _A", " N4 _D"},                     /* C */
      {" N9 _?", " N7 _A", " O6 _A", " N1 _D", " N2 _D", " N3 _A"}, /* G */
      {" N9 _?", " N7 _A", " O6 _A", " N1 _D", " N3 _A"},           /* I */
      {" N1 _?", " O2 _A", " N3 _D", " O4 _A"},                     /* T */
      {" N1 _?", " O2 _A", " N3 _D", " O4 _A"}                      /* U */
  };
  long i, inum = -1, jnum = -1, num = 6;

  if ((pchar = strchr(cmn_base, toupper((int)basei))) != NULL)
    inum = pchar - cmn_base;
  if ((pchar = strchr(cmn_base, toupper((int)basej))) != NULL)
    jnum = pchar - cmn_base;
  if (inum >= 0 && jnum >= 0) { /* both are one of the six bases */
    for (i = 0; i < num; i++) { /* check if backbone atom */
      if (!strncmp(bb_da[i], hb_atom1, 4))
        ia = bb_da[i][5];
      if (!strncmp(bb_da[i], hb_atom2, 4))
        ja = bb_da[i][5];
    }
    if (!ia)
      for (i = 0; i < num; i++)
        if (base_da[inum][i] && !strncmp(base_da[inum][i], hb_atom1, 4)) {
          ia = *(base_da[inum][i] + 5);
          break;
        }
    if (!ja)
      for (i = 0; i < num; i++)
        if (base_da[jnum][i] && !strncmp(base_da[jnum][i], hb_atom2, 4)) {
          ja = *(base_da[jnum][i] + 5);
          break;
        }
    if (ia && ja) {
      sprintf(da, "%c%c", ia, ja);
      if (num_strmatch(da, da_type, 0, 6))
        hbatom_type = '-';
    }
  }

  return hbatom_type;
}

long asym_idx(char *asym, char atoms_list[NELE][3], long dft_lval) {
  long i;

  for (i = 0; i < NELE; i++)
    if (!strcmp(asym, atoms_list[i]))
      return i;

  return dft_lval;
}

void atom_info(long idx, char atoms_list[NELE][3], double *covalence_radii,
               double *vdw_radii) {
  static char *ALIST[NELE] = {UNKATM, " C", " O", " H", " N", " S",
                              " P",   " F", "CL", "BR", " I", "SI"};
  static double CRADII[NELE] = {
      1.666, 0.762, 0.646, 0.352, 0.689, 1.105,
      1.000, 0.619, 1.022, 1.183, 1.378, 1.105 /* XX original 1.200 */
  };
  static double VRADII[NELE] = {2.00, 1.70, 1.52, 1.20, 1.55, 1.80,
                                1.80, 1.47, 1.75, 1.85, 1.98, 2.10};
  long i;

  if (idx == 1)
    for (i = 0; i < NELE; i++)
      strcpy(atoms_list[i], ALIST[i]);
  else if (idx == 2)
    for (i = 0; i < NELE; i++)
      covalence_radii[i] = CRADII[i];
  else if (idx == 3)
    for (i = 0; i < NELE; i++)
      vdw_radii[i] = VRADII[i];
  else
    fatal("wrong options for <atom_info>: should be 1, 2 or 3\n");
}

/* get 2-letter atomic symbol from atom name in PDB file */
void aname2asym(const char *aname0, char *my_asym, long num_sa,
                char **atomlist) {
  char aname[BUF512];
  long i, unknown, NAME_LEN = 4;

  strcpy(aname, aname0);

  /* non-alphabets changed to '.': e.g., " O1P": ".O.P"; " N1 ": ".N.." */
  for (i = 0; i < NAME_LEN; i++)
    if (!isalpha((int)aname[i]))
      aname[i] = '.';

  for (i = 1; i <= num_sa; i++)
    if (str_pmatch(atomlist[i], aname))
      break;

  if (i > num_sa) {
    unknown = is_equal_string(aname, ".UNK");

    if ((aname[0] != '.') && (aname[1] != '.') && (aname[2] == '.') &&
        (aname[3] == '.')) { /* as in 'NA..' */
      my_asym[0] = aname[0];
      my_asym[1] = aname[1];
      my_asym[2] = '\0';
    } else if ((aname[0] == '.') && (aname[1] != '.') && !unknown) {
      my_asym[0] = ' ';
      my_asym[1] = aname[1];
      my_asym[2] = '\0';
    } else if (aname[0] == 'H')
      strcpy(my_asym, " H");
    else
      strcpy(my_asym, UNKATM);

    if (!unknown) {
      fprintf(stderr, "no matching entry for atom name [%s] (%s) in '%s'\n",
              aname0, aname, ATOM_FILE);
      fprintf(stderr, "\tnow it is set as '%s'\n", my_asym);
      fprintf(stderr, "\tcheck and update file $X3DNA/config/atomlist.dat\n");
    }
  } else
    strcpy(my_asym, atomlist[i] + NAME_LEN);
}

/* get atom index for calculating bond linkage */
void atom_idx(long num, char **AtomName, char **Miscs, long *idx) {
  char pdb_asym[3], my_asym[3], atoms_list[NELE][3];
  long i, k, bad_pdb_asym = false;

  atom_info(1, atoms_list, NULL, NULL);

  for (i = 1; i <= num; i++) {
    strcpy(pdb_asym, UNKATM); /* atomic symbol from PDB data file */
    k = false;                /* default to no atomic symbol */

    if (Miscs && strlen(Miscs[i]) >= 27 && !str_pmatch(Miscs[i] + 25, "  ")) {
      strncpy(pdb_asym, Miscs[i] + 25, 2); /* with this info. */
      pdb_asym[2] = '\0';

      if (is_equal_string(pdb_asym, " D"))
        strcpy(pdb_asym, " H");

      if (num_strmatch(pdb_asym, Gvars.ATOM_NAMES, 0, Gvars.NUM_ELE))
        k = true;
      else
        bad_pdb_asym = true; /* for overall checking */
    }

    if (k && !bad_pdb_asym)
      strcpy(my_asym, pdb_asym);
    else /* atomic symbol deduced from atom name */
      aname2asym(AtomName[i], my_asym, Gvars.NUM_SATOM, Gvars.ATOMLIST);

    idx[i] = asym_idx(my_asym, atoms_list, 0);

    /* fprintf(stderr, "%ld\t%s\t%s\t%ld\n", i, AtomName[i], my_asym, idx[i]);
     */
  }

  if (bad_pdb_asym)
    fprintf(stderr, "\tPDB with illegal atomic symbol in columns #77-78\n");
}

/* get linkage information based on inter-atomic distance */
void get_bonds(long num, char **AtomName, double **xyz, long num_residue,
               long *RY, long **seidx, long **connect) {
  long i, nbond, nbond_estimated, *idx, **linkage;

  idx = lvector(1, num);
  atom_idx(num, AtomName, NULL, idx);

  nbond_estimated = lround(NBOND_FNUM * NUM_RESIDUE_ATOMS);
  linkage = lmatrix(1, nbond_estimated, 1, 2);
  for (i = 1; i <= num_residue; i++) /* linkage within each residue */
    if (RY[i] >= 0) {                /* a base residue */
      nbond = 0;
      atom_linkage(seidx[i][1], seidx[i][2], idx, xyz, NULL, NULL,
                   nbond_estimated, &nbond, linkage);
      lkg2connect(AtomName, seidx[i][1], seidx[i][2], nbond, linkage, connect);
    }

  free_lvector(idx, 1, num);
  free_lmatrix(linkage, 1, nbond_estimated, 1, 2);
}

/* get atom linkage using covalent radii criterion. alt. position & chain info.
 */
void atom_linkage(long ib, long ie, long *idx, double **xyz, char **Miscs,
                  char *ChainID, long nbond_estimated, long *nbond,
                  long **linkage) {
  char a1, a2, c1, c2;
  double dst, covalence_radii[NELE];
  long j, k;

  /* Bond criteria: 1.15 * (rA + rB)
     RasMol: 0.56 + (rA + rB)   MacroModel: 1.25 * (rA + rB) */

  atom_info(2, NULL, covalence_radii, NULL);
  for (j = ib; j <= ie - 1; j++) {
    a1 = (Miscs == NULL) ? ' ' : Miscs[j][1]; /* atom alternative position */
    c1 = (ChainID == NULL) ? '_' : ChainID[j];
    for (k = j + 1; k <= ie; k++) {
      a2 = (Miscs == NULL) ? ' ' : Miscs[k][1]; /* atom alternative position */
      c2 = (ChainID == NULL) ? '_' : ChainID[k];
      if (a1 != ' ' && a2 != ' ' && a1 != a2) /* same model: alt. position */
        continue;
      if (c1 != c2) /* same chain ID */
        continue;
      dst = BOND_FACTOR * (covalence_radii[idx[j]] + covalence_radii[idx[k]]);
      if (within_limits(xyz[k], xyz[j], 0, dst)) {
        if (++*nbond > nbond_estimated)
          fatal("too many linkages\n");
        else {
          linkage[*nbond][1] = j;
          linkage[*nbond][2] = k;
        }
      }
    }
  }
}

/* from bond linkage to PDB connection table: SLOW for large structure rr0033 */
void lkg2connect(char **AtomName, long ib, long ie, long nbond, long **linkage,
                 long **connect) {
  long i, ilink, j;
  long idx[7];

  for (i = ib; i <= ie; i++) {
    ilink = 0;
    for (j = 1; j <= nbond; j++) {
      if (i == linkage[j][1]) {
        if (++ilink > 6) {
          fprintf(stderr, "atom <%ld: [%s]> has over 6 bonds\n", i,
                  AtomName[i]);
          break;
        } else
          idx[ilink] = linkage[j][2];
      }
      if (i == linkage[j][2]) {
        if (++ilink > 6) {
          fprintf(stderr, "atom <%ld: [%s]> has over 6 bonds\n", i,
                  AtomName[i]);
          break;
        } else
          idx[ilink] = linkage[j][1];
      }
    }
    if (ilink > 6)
      ilink = 6; /* maximum six bonds */
    if (ilink > 0) {
      for (j = 1; j <= ilink; j++)
        connect[i][j] = idx[j];
      connect[i][7] = ilink;
      lsort(ilink, connect[i], idx); /* ignore idx */
    }
  }
}

/* initialize matrix htm_water as follows:
 * row#1: c0 for # of atoms; c1-num for sequential residue # in the PDB file
 * row#2: c0 for # of residues; c1-num for connected base sequential residue #
 * row#3: c0 for hydration indication; c1-num for atom index
 * row#4: c0 for # of waters; c1--num_H2O for corresponding residue # */
void init_htm_water(long waters, long num, long num_residue, long *idx,
                    long **htm_water) {
  long i;

  htm_water[1][0] = num;         /* total # of atoms */
  htm_water[2][0] = num_residue; /* total # of residue */
  htm_water[3][0] = waters;      /* if to check for hydration */

  for (i = 1; i <= num; i++) /* make a copy of atom index */
    htm_water[3][i] = idx[i];
}

/* make each HETATM record, including water, to its connected base residue */
void identify_htw(long num_residue, long **seidx, long *RY, char **AtomName,
                  char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                  double **xyz, long **htm_water) {
  static char *WATER[] = {WATER_LIST};
  char a1, a2, c1, c2;
  double dst, covalence_radii[NELE];
  long num_wat = sizeof WATER / sizeof WATER[0] - 1;
  long i, ib, ie, id, j, k = 0, m, num_H2O = 0;

  /* identify ligands and water molecules */
  for (i = 1; i <= num_residue; i++) {
    if (RY[i] >= 0) /* skip base residues */
      continue;
    ib = seidx[i][1];
    ie = seidx[i][2];
    id = ie - ib + 1;
    for (j = ib; j <= ie; j++) {
      if (Miscs[j][0] == 'H') { /* HETATM record */
        htm_water[1][j] = i;
        htm_water[2][j] = -8; /* general HETATM record indicator */
      }
    }

    if (num_strmatch(ResName[ib], WATER, 0, num_wat) /* matching H2O ResName */
        ||
        (id == 1 &&               /* ONLY ONE atom in this residue */
         ((htm_water[3][ib] == 2) /* added nov-03-2006: 2 is Oxygen */
          || !strcmp(AtomName[ib], " O  ") || !strcmp(AtomName[ib], " OW ")))) {
      for (j = ib; j <= ie; j++) { /* possibly with Hs, e.g. NMR */
        htm_water[1][j] = i;
        htm_water[2][j] = -1; /* for water molecule */
      }
      num_H2O++;
      htm_water[4][num_H2O] = i;
    }
  }
  htm_water[4][0] = num_H2O;

  /* This section checks for covalent connection of a non-base residue,
   * e.g., drug molecule, to a base residue. They have to be on the same
   * alternative position and chain ID. Pure H2Os do not count here. */
  atom_info(2, NULL, covalence_radii, NULL); /* cf. <atom_linkage> */
  for (i = 1; i <= num_residue; i++) {
    if (RY[i] >= 0) /* skip base residues */
      continue;
    id = 0; /* check if residue "i" is connected to a base residue */
    for (j = seidx[i][1]; j <= seidx[i][2]; j++) {
      /* assuming NO more than 1 connection to base; only HETATM/H2O */
      if (htm_water[2][j] >= 0) /* -8 or -1 */
        break;
      a1 = Miscs[j][1];
      c1 = ChainID[j];
      for (k = 1; k <= num_residue; k++) {
        if (RY[k] < 0) /* check only for connections to base residue */
          continue;
        for (m = seidx[k][1]; m <= seidx[k][2]; m++) {
          a2 = Miscs[m][1];
          c2 = ChainID[m];
          if (c1 != c2) /* same chain ID */
            break;      /* restricted to residue m */
          if (a1 != ' ' && a2 != ' ' && a1 != a2)
            continue; /* same model: alt. position */
          dst = BOND_FACTOR * (covalence_radii[htm_water[3][j]] +
                               covalence_radii[htm_water[3][m]]);
          if (within_limits(xyz[m], xyz[j], 0, dst)) {
            id = 1;
            goto ID_CNCT;
          }
        }
      }
    }

  ID_CNCT:
    if (id) { /* residue i connected to base residue k */
      for (j = seidx[i][1]; j <= seidx[i][2]; j++)
        htm_water[2][j] = k;
      ib = seidx[i][1]; /* HETATM residue */
      ie = seidx[k][1]; /* contacted base residue */

      /* no longer check for match of residue #: ResSeq[ib] != ResSeq[ie]: 4tna
       */
      if (Miscs[ib][2] !=
          Miscs[ie][2]) /* e.g., 1fjb: 17(697,23) and 17(507,17) */
        fprintf(stderr,
                "%ld(%ld,%ld) and %ld(%ld,%ld) have different"
                " insertion code!\n",
                ResSeq[ib], ib, i, ResSeq[ie], ie, k);
    }
  }
}

/* attach covalently connected residues as in +A etc */
long attached_residues(long inum_base, long *ivec, long *ivec2, long **seidx,
                       double **xyz, long **htm_water, miscPars *misc_pars) {
  long i, iw, ir, j, k, m, n, tnum_res = inum_base;
  long num_residue = htm_water[2][0], num_H2O = htm_water[4][0];

  /* check if there is a covalent bond linkage between hetero
   * group (e.g., a drug molecule) with the base residue. */
  for (i = 1; i <= inum_base; i++) {
    ivec2[i] = ivec[i];

    if (Gvars.ATTACH_RESIDUE ==
        false) /* Per Pascal's request: no metal, or HETATM */
      continue;

    for (j = 1; j <= num_residue; j++) {
      k = seidx[j][1];
      if (htm_water[2][k] == ivec[i])
        ivec2[++tnum_res] = htm_water[1][k];
    }
  }

  if (!htm_water[3][0])
    return tnum_res;

  inum_base = tnum_res;
  n = inum_base + 1; /* starting index for hydration residues */
  for (ir = 1; ir <= num_H2O; ir++) { /* time-consuming part: best way! */
    iw = htm_water[4][ir];            /* hydration index */
    for (m = seidx[iw][1]; m <= seidx[iw][2]; m++) {
      if (htm_water[3][m] != 2) /* not O */
        continue;
      for (i = 1; i <= inum_base; i++) {
        k = ivec2[i];
        for (j = seidx[k][1]; j <= seidx[k][2]; j++) {
          if (!lval_in_set(htm_water[3][j], 1, misc_pars->water_idx[0],
                           misc_pars->water_idx))
            continue;
          if (within_limits(xyz[m], xyz[j], misc_pars->hb_lower,
                            misc_pars->water_dist) &&
              !lval_in_set(iw, n, tnum_res, ivec2)) {
            ivec2[++tnum_res] = iw;
          }
        }
      }
    }
  }

  if (tnum_res > inum_base + 1) { /* sort the hydration residues into order */
    long *idx;

    k = tnum_res - inum_base;
    idx = lvector(1, k);
    lsort(k, ivec2 + inum_base, idx);

    free_lvector(idx, 1, k);
  }

  return tnum_res;
}

/* print out base-pairing information */
void print_pairinfo(long i, long j, char basei, char basej, double *rtn_val,
                    double *chi, miscPars *misc_pars, long **seidx, long *idx,
                    char **AtomName, double **xyz, char *bseq, long detailed,
                    FILE *fp) {
  char antip, bptype[4], hb_info[BUF512];
  long k;

  antip = (rtn_val[35] < 0.0) ? '-' : '+';
  sprintf(bptype, "%c%c%c", bseq[i], antip, bseq[j]);
  fprintf(fp, "             %s  ", bptype);
  for (k = 27; k <= 32; k++) /* bp parameters */
    fprintf(fp, "%8.2f", rtn_val[k]);
  fprintf(fp, "\n");
  get_hbond_ij(i, j, basei, basej, misc_pars, seidx, idx, AtomName, xyz,
               hb_info);
  fprintf(fp, "             %s %s\n", bptype, hb_info);

  if (!detailed)
    return;

  fprintf(fp, "                    %s %s %s\n",
          rtn_val[35] < 0.0 ? "anti-parallel" : "parallel", /* z-axis */
          rtn_val[33] < 0.0 ? "trans" : "cis",              /* x-axis */
          rtn_val[36] < 0.0 ? "trans" : "cis");             /* C1-N vector */
  fprintf(fp, "            ");
  for (k = 33; k <= 35; k++)
    fprintf(fp, "%c", rtn_val[k] > 0.0 ? '+' : '-');
  for (k = 33; k <= 36; k++)
    fprintf(fp, "%9.1f", dot2ang(rtn_val[k]));
  fprintf(fp, "%9.1f%9.1f\n", chi[i], chi[j]);

  fprintf(fp, "     ");
  for (k = 1; k <= 8; k++)
    fprintf(fp, "%8.2f", rtn_val[k]);
  fprintf(fp, "\n");
}

/* once a pair is found, calculate more parameters */
static void calculate_more_bppars(long i, long j, double dir_x, double dir_y,
                                  double dir_z, double **orien, double **org,
                                  char *bseq, double **NC1xyz, double *rtn_val,
                                  long *bpid) {
  char bpi[3];
  long k, l, koffset;
  double zave[4], dNN_vec[4], pars[7], **r1, **r2, **mst;

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mst = dmatrix(1, 3, 1, 3);

  for (k = 1; k <= 9; k++) {
    rtn_val[k + 8] = orien[i][k];  /* base I xyz axes */
    rtn_val[k + 17] = orien[j][k]; /* base II xyz axes */
  }
  for (k = 1; k <= 3; k++) {
    koffset = (k - 1) * 3;
    for (l = 1; l <= 3; l++) {
      r1[l][k] = orien[i][koffset + l];
      r2[l][k] = (k == 1 || dir_z > 0)
                     ? /* keep x, reverse y & z if anti-parallel */
                     orien[j][koffset + l]
                     : -orien[j][koffset + l];
    }
  }
  bpstep_par(r2, org[j], r1, org[i], pars, mst, &rtn_val[5]);
  for (k = 1; k <= 6; k++) /* bp parameters in columns 27-32 */
    rtn_val[26 + k] = pars[k];

  sprintf(bpi, "%c%c", toupper((int)bseq[i]), toupper((int)bseq[j]));
  *bpid = -1; /* assumed to be non-WC */
  if (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0) {
    check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
    if (*bpid == 2)
      rtn_val[5] -= 2.0; /* WAS: -1.5 (bonus for WC pair) */
  }

  rtn_val[33] = dir_x;
  rtn_val[34] = dir_y;
  rtn_val[35] = dir_z;

  if (NC1xyz[i][7] > 0 && NC1xyz[j][7] > 0) {
    ddxyz(NC1xyz[i], NC1xyz[i] + 3, zave);
    ddxyz(NC1xyz[j], NC1xyz[j] + 3, dNN_vec);
    vec_norm(zave);
    vec_norm(dNN_vec);
    rtn_val[36] = dot(zave, dNN_vec);
  } else
    rtn_val[36] = EMPTY_NUMBER;

  free_dmatrix(r1, 1, 3, 1, 3);
  free_dmatrix(r2, 1, 3, 1, 3);
  free_dmatrix(mst, 1, 3, 1, 3);
}

/* for judging the quality of a pair. mostly used to find the "best" pair */
static double adjust_pairQuality(long i, long j, char basei, char basej,
                                 long **seidx, long *idx, char **AtomName,
                                 double **xyz, miscPars *misc_pars) {
  double dval;
  long k, num_hb, num_good_hb = 0, **num_list;

  num_list = lmatrix(1, BUF512, 0, 3);
  hb_numlist(i, j, basei, basej, seidx, idx, AtomName, xyz, misc_pars, &num_hb,
             num_list);
  for (k = 1; k <= num_hb; k++) {
    if (num_list[k][0]) /* not canonical D-A type */
      continue;
    dval = num_list[k][3] / MFACTOR;
    if (dval_in_range(dval, 2.5, 3.5)) /* assumed good H-bonding distance */
      num_good_hb++;
  }
  free_lmatrix(num_list, 1, BUF512, 0, 3);

  if (num_good_hb >= 2) /* make it works for i54:55 vs j238:239 of 1zla */
    return -3.0;
  else
    return -num_good_hb;
}

/* Checking if two bases form a pair according to several criteria
 * rtn_val[RTNNUM]:
 *       d, dv, angle, dNN, dsum, bp-org,  x1,  y1,   z1,   x2,   y2,   z2
 * col#  1   2    3     4     5    6-8    9-11 12-14 15-17 18-20 21-23 24-26
 * bpid: 0: not-paired; +1: WC geometry; +2: WC pair; -1: other cases
 *       bp-pars,   bp_relative_orientation, C1_N_relative_orientation
 * col#   27-32        33-35 (dx, dy, dz)              36
 *
 * There are also possibilities when the pair is mostly maintained by
 * backbone + backbone; backbone + base H-bonds. e.g. A512 + C637
 * in rr0033 (chain 0). No direct base-to-base H-bonds ************** */
void check_pair(long i, long j, char *bseq, long **seidx, double **xyz,
                double **NC1xyz, double **orien, double **org, long *idx,
                char **AtomName, miscPars *misc_pars, double *rtn_val,
                long *bpid, long **ring_atom, long network) {
  double dir_x, dir_y, dir_z;
  double dorg[4], oave[4], zave[4], dNN_vec[4];
  long cdns, m, n, num_base_hb = 0, num_o2_hb = 0;

  *bpid = 0; /* default as not-paired */
  if (i == j)
    return; /* same residue */

  get_bp_zoave(i, j, orien, org, oave, zave);

  ddxyz(org[i], org[j], dorg);
  ddxyz(NC1xyz[i], NC1xyz[j], dNN_vec);

  rtn_val[1] = veclen(dorg);               /* distance between origins */
  dir_x = dot(&orien[i][0], &orien[j][0]); /* relative x direction */
  dir_y = dot(&orien[i][3], &orien[j][3]); /* relative y direction */
  dir_z = dot(&orien[i][6], &orien[j][6]); /* relative z direction */

  rtn_val[2] = fabs(dot(dorg, zave)); /* dv: projection onto mean normal */
  rtn_val[3] = z1_z2_angle_in_0_to_90(
      &orien[i][6], &orien[j][6]); /* angle between base normals */
  rtn_val[4] = veclen(dNN_vec);    /* RN9-YN1 distance */
  rtn_val[5] = rtn_val[1] + 2.0 * rtn_val[2] + rtn_val[3] / 20.0;
  /* WAS 25. changed based on 1id3 i25:26/j267:268; i107:8/j185:6 */

  if (network) { /* check if two bases in pairing network */
#if 0
        fprintf(stderr, "i: %ld\tj: %ld\tangle: %g\n", i, j, rtn_val[3]);
#endif
    if (dval_in_range(rtn_val[3], misc_pars->min_plane_angle,
                      misc_pars->max_plane_angle) &&
        dval_in_range(rtn_val[4], misc_pars->min_dNN, misc_pars->max_dNN)
        /* -- more strict than for best-pair below; otherwise,
           np_recipes/R5_pentaplets will have different result at least for
           'multiplets.pdb'
                    && get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVERLAP)
        */
        && get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVERLAP &&
        (get_oarea(i, j, ring_atom, org[i], &orien[i][6], xyz, 0) < OVERLAP) &&
        (get_oarea(i, j, ring_atom, org[j], &orien[j][6], xyz, 0) < OVERLAP))
      *bpid = -1; /* WAS 1, Oct. 30, 2006 */
    rtn_val[5] += adjust_pairQuality(i, j, bseq[i], bseq[j], seidx, idx,
                                     AtomName, xyz, misc_pars);
    return;
  }

  cdns = (dval_in_range(rtn_val[1], misc_pars->min_dorg, misc_pars->max_dorg) &&
          dval_in_range(rtn_val[2], misc_pars->min_dv, misc_pars->max_dv) &&
          dval_in_range(rtn_val[3], misc_pars->min_plane_angle,
                        misc_pars->max_plane_angle) &&
          dval_in_range(rtn_val[4], misc_pars->min_dNN, misc_pars->max_dNN));
  /* no-overlap with reference to middle, lower (i) or upper (j) frame
      if (cdns &&
          (get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVERLAP) &&
          (get_oarea(i, j, ring_atom, org[i], &orien[i][6], xyz, 0) < OVERLAP)
     && (get_oarea(i, j, ring_atom, org[j], &orien[j][6], xyz, 0) < OVERLAP)) {
  */
  if (cdns && get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVERLAP) {
    for (m = seidx[i][1]; m <= seidx[i][2]; m++)
      for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
        if (!within_limits(xyz[m], xyz[n], misc_pars->hb_lower,
                           misc_pars->hb_dist1))
          continue;

        if (is_baseatom(AtomName[m]) && is_baseatom(AtomName[n]) &&
            good_hbatoms(misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]))
          num_base_hb++;

        if (is_equal_string(AtomName[m], " O2'") ||
            is_equal_string(AtomName[n], " O2'"))
          num_o2_hb++;
      }

    if ((misc_pars->min_base_hb && (num_base_hb >= misc_pars->min_base_hb)) ||
        (!misc_pars->min_base_hb && (num_o2_hb || num_base_hb))) {
      calculate_more_bppars(i, j, dir_x, dir_y, dir_z, orien, org, bseq, NC1xyz,
                            rtn_val, bpid);
      rtn_val[5] += adjust_pairQuality(i, j, bseq[i], bseq[j], seidx, idx,
                                       AtomName, xyz, misc_pars);
    }
  }
}

/* get the xyz coordinates of O3' and P atoms */
void o3_p_xyz(long ib, long ie, char *aname, char **AtomName, double **xyz,
              double *o3_or_p, long idx) {
  long i;

  i = find_1st_atom(aname, AtomName, ib, ie, "");
  if (i) { /* with O3'/P atom */
    cpxyz(xyz[i], o3_or_p + idx - 4);
    o3_or_p[idx] = 1.0;
  } else /* without O3'/P atom */
    o3_or_p[idx] = -1.0;
}

/* check if atomname is a nucleotide base atom: ASSUMING NA residue */
long is_baseatom(char *atomname) {
  if (is_equal_string(atomname, " C5M")) /* C5M of T */
    return true;

  if (atomname[0] == ' ' && strchr("HP", atomname[1]) == NULL /* like " N1 " */
      && isdigit((int)atomname[2]) && atomname[3] == ' ')
    return true;

  return false;
}

static long glyco_N(long isR, long ib, long ie, char b, char **AtomName,
                    char **ResName, long C1prime, double **xyz) {
  char *a;
  double d, dm = 9999.0, *c1xyz;
  long k, km, num = 0;

  if (isR) { /* purine */
    k = find_1st_atom(" N9 ", AtomName, ib, ie, "");
    if (k)
      return k;
  } else {                 /* pyrimidine */
    if (strchr("Pp", b)) { /* pseudo-uridine */
      k = find_1st_atom(" C5 ", AtomName, ib, ie, "");
      if (k)
        return k;
    } else { /* regular pyrimidine */
      k = find_1st_atom(" N1 ", AtomName, ib, ie, "");
      if (k)
        return k;
    }
  }

  assert(!k);
  fprintf(stderr, "Cannot identify RN9/YN1 in [%s] -- ", ResName[ib]);

  if (C1prime) { /* find the shorest distance */
    c1xyz = xyz[C1prime];
    km = false;
    for (k = ib; k <= ie; k++) {
      a = AtomName[k];
      if (!is_baseatom(a))
        continue;
      d = p1p2_dist(c1xyz, xyz[k]);
      if (d < dm) {
        dm = d;
        km = k;
      }
    }
  }

  if (dm <= BOND_UPPER_LIMIT) {
    fprintf(stderr, "use atom [%s] instead\n", AtomName[km]);
    return km;
  }

  /* last try! */
  km = false;
  for (k = ib; k <= ie; k++) {
    a = AtomName[k];
    if ((isR && strchr(a, '9')) || (!isR && strchr(a, '1'))) {
      num++;
      km = k;
    }
  }

  if (num == 1) {
    fprintf(stderr, "[i] using atom [%s] instead\n", AtomName[km]);
    return km;
  }

  fatal("stop!\n");

  return 0;
}

/* get base information for locating possible pairs later */
void base_info(long num_residue, char *bseq, long **seidx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, double **orien, double **org,
               double **NC1xyz, double **o3_p) {
  char BDIR[BUF1K];
  long i, ib, ie, C1prime, N;

  get_BDIR(BDIR, "Atomic_A.pdb");

  /* get the reference frame for each base */
  base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq,
             Miscs, xyz, BDIR, orien, org);

  for (i = 1; i <= num_residue; i++) {
    ib = seidx[i][1];
    ie = seidx[i][2];
    if (RY[i] >= 0) { /* a base */
      C1prime = find_1st_atom(" C1'", AtomName, ib, ie, "");
      if (C1prime) {
        NC1xyz[i][7] = 1.0;                 /* C1' atom exists */
        cpxyz(xyz[C1prime], NC1xyz[i] + 3); /* 4-6 */
      }
      N = glyco_N(RY[i], ib, ie, bseq[i], AtomName, ResName, C1prime, xyz);
      cpxyz(xyz[N], NC1xyz[i]); /* RN9/YN1 surely exist */
      o3_p_xyz(ib, ie, " O3'", AtomName, xyz, o3_p[i], 4);
      o3_p_xyz(ib, ie, " P  ", AtomName, xyz, o3_p[i], 8);
    }
  }
}

void help3dna_usage(char *program_name) {
  help3dna(program_name);
  contact_msg(0);
}

/* Easy 3DNA help message read from a text file */
void help3dna(char *program_name) {
  char BDIR[BUF1K], str[BUF512], pname[BUF512], *prefix = "        ";
  long nlen, nchar = 75, found_help = 0;
  FILE *fp;

  get_BDIR(BDIR, HELP3DNA);
  strcat(BDIR, HELP3DNA);
  fp = open_file(BDIR, "r");

  print_sep(stderr, '=', nchar);
  strcpy(pname, program_name); /* program_name is a constant */
  nlen = upperstr(pname);      /* to upper case */
  while (fgets(str, sizeof str, fp) != NULL) {
    if (str[0] != '<')
      continue;
    upperstr(str);
    if (strncmp(str + 1, pname, nlen))
      continue; /* not match */
    found_help = 1;
    while (fgets(str, sizeof str, fp) != NULL) {
      if (str[0] == '<') {
        upperstr(str);
        if (str[1] != '/' || (strncmp(str + 2, pname, nlen)))
          fatal("error in help format: <tag> ... </tag>\n");
        found_help = 2;
        goto FINISHED;
      } else {
        if (str[0] != '#')            /* not-comment */
          fprintf(stderr, "%s", str); /* str includes \n */
      }
    }
  }

FINISHED:
  if (!found_help)
    fprintf(stderr, "No help found for program <%s>\n", program_name);
  else if (found_help != 2)
    fprintf(stderr, "Warning: no matching tag found\n");
  else {
    fprintf(stderr, "AUTHOR\n%s%s\n\n", prefix, Gvars.X3DNA_VER);
    fprintf(stderr, "Please post questions/comments on the 3DNA Forum:"
                    " http://forum.x3dna.org/\n");
    fprintf(stderr, "Please check 'http://x3dna.org/citations' on how to cite"
                    " 3DNA --- THANKS!\n");
  }

  print_sep(stderr, '=', nchar);

  close_file(fp);
}

/* Delete all H-atoms and write the rest coordinates in PDB format */
void delH_pdbfile(char *inpfile, char *outfile) {
  char *ChainID, *p, **AtomName, **ResName, **Miscs;
  double **xyz;
  long i, k = 0, num, *ResSeq, *idx;
  FILE *fp;

  /* read in the PDB file */
  num = number_of_atoms(inpfile, 1, "*");
  AtomName = cmatrix(1, num, 0, 4);
  ResName = cmatrix(1, num, 0, 3);
  ChainID = cvector(1, num);
  ResSeq = lvector(1, num);
  xyz = dmatrix(1, num, 1, 3);
  Miscs = cmatrix(1, num, 0, NMISC);
  read_pdb(inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1,
           "*");

  idx = lvector(1, num);
  atom_idx(num, AtomName, Miscs, idx); /* atom index to mark H */

  fp = open_file(outfile, "w");
  print_pdb_title(inpfile, "*", fp);
  for (i = 1; i <= num; i++) {
    if (idx[i] == 3) /* H-atom */
      continue;
    p = Miscs[i];
    fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
            (p[0] == 'A') ? "ATOM  " : "HETATM", ++k, AtomName[i], p[1],
            ResName[i], ChainID[i], ResSeq[i], p[2], xyz[i][1], xyz[i][2],
            xyz[i][3], p + 3);
  }
  fprintf(fp, "END\n");
  close_file(fp);

  free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
  free_lvector(idx, 1, num);
}

void contact_msg(long prt_msg) {
  if (prt_msg) {
    fprintf(stderr, "%s\n", Gvars.X3DNA_VER);
    help3dna("contact_info");
  }
  exit(1);
}

/* partial command-line option match from the beginning */
long str_pmatch(char *str, char *sstr) {
  return !strncmp(str, sstr, strlen(sstr));
}

/* case-insensitive command-line option match from the beginning */
long case_str_pmatch(char *str, char *sstr) {
  return !case_strncmp(str, sstr, strlen(sstr));
}

static long isEMPTY_or_isNAN_or_isINF(char *str) {
  lowerstr(str); /* to lower case for consistency: NAN/nan etc */

  /* empty field/missing value etc: \t\t; NaN; Inf/-Inf */
  if (*str == '\0' || strstr(str, "nan") || strstr(str, "inf"))
    return true;
  else
    return false;
}

/* convert string "str" to a double value, with error checking */
double cvt2double(char *str) {
  char *endp, *p = trim(str);
  double d;

  if (isEMPTY_or_isNAN_or_isINF(p))
    return XBIG; /* as an impossible value */

  errno = 0;
  d = strtod(p, &endp);
  if (*endp != '\0' || errno == ERANGE)
    return XBIG; /* for NA etc */
  else
    return d;
}

/* convert string "str" to a long value, with error checking */
long cvt2long(char *str) {
  char *endp, *p = trim(str);
  long d;

  if (isEMPTY_or_isNAN_or_isINF(p))
    return LONG_MAX; /* as an impossible value */

  errno = 0;
  d = strtol(p, &endp, 10);
  if (*endp != '\0' || errno == ERANGE)
    return LONG_MAX; /* for NA etc */
  else
    return d;
}

long equalsign_pos(char *str) {
  char *ptr;

  ptr = strchr(str, '='); /* check from the beginning */
  if (ptr == NULL)
    fatal("wrong format for var=value pair [%s]\n", str);
  return ptr - str + 1;
}

/* extract a long numerical value from command line */
long get_lvalue(char *str, long vmin, long vmax) {
  long npos, val;

  npos = equalsign_pos(str);

  val = cvt2long(str + npos);
  if (val == LONG_MAX)
    fatal("wrong option [%s]: not a valid integer value\n", str);

  if (!lval_in_range(val, vmin, vmax))
    fatal("invalid option [%s]: value %ld out of range [%ld %ld]\n", str, val,
          vmin, vmax);

  return val;
}

/* extract a double numerical value from command line */
double get_dvalue(char *str, double vmin, double vmax) {
  long npos;
  double val;

  npos = equalsign_pos(str);

  val = cvt2double(str + npos);
  if (val > XBIG_CUTOFF)
    fatal("wrong option [%s]: not a valid numerical value\n", str);

  if (!dval_in_range(val, vmin, vmax))
    fatal("invalid option [%s]: value %f out of range [%f %f]\n", str, val,
          vmin, vmax);

  return val;
}

/* extract a string from command line, substitute ~ to $HOME */
void get_strvalue(char *str, char *dst, long expand_tilde) {
  char *p;
  long npos;

  if (strlen(str) > BUF512)
    fatal("command line option too long: %.36s...\n", str);

  npos = equalsign_pos(str);
  if (expand_tilde && str[npos] == '~') {
    p = getenv("HOME");
    if (p == NULL) {
      fprintf(stderr, "no environment variable HOME defined!\n");
      strcpy(dst, str + npos);
    } else
      sprintf(dst, "%s%s", p, str + npos + 1);
  } else
    strcpy(dst, str + npos);
}

/* convert all occurrences of character 'c1' to 'c2' in string 'str' */
void cvtstr_c1toc2(char *str, char c1, char c2) {
  char *p = str;

  while (*p) {
    if (*p == c1)
      *p = c2;
    p++;
  }
}

/* assuming z1 and z2 are already normalized */
double z1_z2_angle_in_0_to_90(double *z1, double *z2) {
  double dircos = dot(z1, z2);

  return 90.0 - fabs(dot2ang(dircos) - 90.0); /* absolute value */
}

void do_nothing(void) { return; }

void skip_lines(long num, FILE *fp) {
  char *p0;
  long i;

  for (i = 1; i <= num; i++) {
    p0 = my_getline(fp);
    if (p0 == NULL)
      fatal("no <%ld> lines to skip!\n", num);
    free(p0);
  }
}

void print_frame(FILE *fp, double *O, double **R) {
  fprintf(fp, "%10.4f %10.4f %10.4f  # origin\n", O[1], O[2], O[3]);
  fprintf(fp, "%10.4f %10.4f %10.4f  # x-axis\n", R[1][1], R[2][1], R[3][1]);
  fprintf(fp, "%10.4f %10.4f %10.4f  # y-axis\n", R[1][2], R[2][2], R[3][2]);
  fprintf(fp, "%10.4f %10.4f %10.4f  # z-axis\n", R[1][3], R[2][3], R[3][3]);
}

/// @brief - analyze.cpp

using namespace std;

#define SIMPLE_BP_LONG_AXIS_RN9_YN1 2
#define SIMPLE_STEP_HELICAL_PARS 4

typedef struct {
  char torsion[BUF512];
  long istart;
  long istep;
  long icnt;
  long waters;
  long bz;          /* for B-Z junction */
  long ring;        /* for base ring center and normal vector */
  long simple_pars; /* simplified base-pair (YC6--RC8) and step (C1'--C1')
                       parameters */
  long abi;         /* Stephen Harvey's A/B index */
  long circular;
} struct_args;

static void get_ring_center_normal(long ds, long num_bp, long **pair_num,
                                   char **bp_seq, long **seidx, long *RY,
                                   char **AtomName, char **ResName,
                                   char *ChainID, long *ResSeq, char **Miscs,
                                   double **xyz, FILE *fp) {
  static char *RingAtom[] = {RA_LIST};
  char BDIR[BUF512], idmsg[BUF512], sidmsg[BUF512], spdb[BUF512];
  char *sChainID, **sAtomName, **sResName, **sMiscs, *fmt = " %9.3f";

  double orgi[4], cxyz[4], nvec[4], rmsd;
  double **eRing_xyz, **fitted_xyz, **sRing_xyz, **sxyz, **R;

  long i, ib, ie, j, k, rnum, RingAtom_num, RA_NUM;
  long exp_katom, nmatch, snum, std_katom;
  long *sResSeq;

  get_BDIR(BDIR, "Atomic_A.pdb");
  RA_NUM = sizeof RingAtom / sizeof RingAtom[0];

  eRing_xyz = dmatrix(1, RA_NUM, 1, 3);
  sRing_xyz = dmatrix(1, RA_NUM, 1, 3);
  sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
  sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
  sChainID = cvector(1, NUM_RESIDUE_ATOMS);
  sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
  sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
  sMiscs = cmatrix(1, NUM_RESIDUE_ATOMS, 0, NMISC);
  fitted_xyz = dmatrix(1, RA_NUM, 1, 3);
  R = dmatrix(1, 3, 1, 3);

  for (i = 1; i <= ds; i++) {
    print_sep(fp, '*', 76);
    if (ds == 2)
      fprintf(fp,
              "Strand %s base ring center (Ox, Oy, Oz) and normal vector"
              " (Nx, Ny, Nz) in\n   the coordinate system"
              " of the given structure\n\n",
              (i == 1) ? "I" : "II");
    else
      fprintf(fp,
              "Base ring center (Ox, Oy, Oz) and normal vector (Nx, Ny, Nz) in"
              "\n   the coordinate system of the given structure\n\n");

    fprintf(fp, "    base        Ox        Oy        Oz        Nx        Ny    "
                "    Nz\n");

    for (j = 1; j <= num_bp; j++) {
      rnum = pair_num[i][j];
      ib = seidx[rnum][1];
      ie = seidx[rnum][2];
      get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);

      RingAtom_num = (RY[rnum] == 1) ? RA_NUM : RA_NUM - 3;
      set_std_base_pdb(BDIR, false, bp_seq[i][j], spdb);
      snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                      sMiscs, 1, "*");
      sprintf(sidmsg, "in standard base: %s", spdb);

      nmatch = 0;
      for (k = 0; k < RingAtom_num; k++) {
        exp_katom = find_1st_atom(RingAtom[k], AtomName, ib, ie, idmsg);
        std_katom = find_1st_atom(RingAtom[k], sAtomName, 1, snum, sidmsg);
        if (exp_katom && std_katom) {
          ++nmatch;
          cpxyz(xyz[exp_katom], eRing_xyz[nmatch]);
          cpxyz(sxyz[std_katom], sRing_xyz[nmatch]);
        }
      }
      rmsd = ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, orgi);
      UNUSED_PARAMETER(rmsd);

      ave_dmatrix(eRing_xyz, nmatch, 3, cxyz);
      for (k = 1; k <= 3; k++) /* column-wise */
        nvec[k] = R[k][3];

      fprintf(fp, " %4ld %c   ", j, bp_seq[i][j]);
      for (k = 1; k <= 3; k++)
        fprintf(fp, fmt, cxyz[k]); /* center */
      for (k = 1; k <= 3; k++)
        fprintf(fp, fmt, nvec[k]); /* base normal */
      fprintf(fp, "\n");
    }
  }

  free_dmatrix(eRing_xyz, 1, RA_NUM, 1, 3);
  free_dmatrix(sRing_xyz, 1, RA_NUM, 1, 3);
  free_cmatrix(sAtomName, 1, NUM_RESIDUE_ATOMS, 0, 4);
  free_cmatrix(sResName, 1, NUM_RESIDUE_ATOMS, 0, 3);
  free_cvector(sChainID, 1, NUM_RESIDUE_ATOMS);
  free_lvector(sResSeq, 1, NUM_RESIDUE_ATOMS);
  free_dmatrix(sxyz, 1, NUM_RESIDUE_ATOMS, 1, 3);
  free_cmatrix(sMiscs, 1, NUM_RESIDUE_ATOMS, 0, NMISC);
  free_dmatrix(fitted_xyz, 1, RA_NUM, 1, 3);
  free_dmatrix(R, 1, 3, 1, 3);
}

static void simple_bp_pars(double **orien, double **org, long NNy, long **c6_c8,
                           long **chi, long num_bp, char **bp_seq,
                           long *WC_info, double **xyz, FILE *fp) {
  char wc, str[BUF512], tmp[BUF32];
  char *bstr = "      ----", *fmt = "%10.2f";
  long i, j, k, Ly, Ry, num = 0;
  double buckle, propeller, norm, r, angle;
  double dx[4], dy[4], dz[4], pars[7];
  double x1[4], y1[4], z1[4], x2[4], y2[4], z2[4];
  double morg[4], o1[4], o2[4], dd[4];
  double **r1, **r2, **mst, **parcln;

  print_sep(fp, '-', 76);
  fprintf(fp, "Simple base-pair parameters based on %s vectors\n",
          NNy ? "RN9--YN1" : "RC8--YC6");
  fprintf(fp, "      bp        Shear    Stretch   Stagger    Buckle  Propeller "
              " Opening\n");

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mst = dmatrix(1, 3, 1, 3);
  parcln = dmatrix(1, num_bp, 1, 6);

  for (i = 1; i <= num_bp; i++) {
    /* Here: r1 refers to the base on strand II; r2 for strand I base.
     * For NEGATIVE orientaion, the r1 y/z-axes are already flipped */
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    bpstep_par(r1, o1, r2, o2, pars, mst, morg);
    mtx_2_x_y_z(r1, x1, y1, z1);
    mtx_2_x_y_z(r2, x2, y2, z2);
    mtx_2_x_y_z(mst, dx, dy, dz); /* keep dz as is */

    angle = z1_z2_angle_in_0_to_90(z1, z2);
    wc = (WC_info[i] == 2) ? ' ' : '*';
    sprintf(str, "%c %4ld %c%c%c ", wc, i, bp_seq[1][i], bp_seq[0][i],
            bp_seq[2][i]);

    if (NNy) { /* RN9--YN1 */
      k = (i - 1) * 4;
      Ly = chi[1][k + 3];
      Ry = chi[2][k + 3];
    } else { /* RC8--YC6 */
      Ly = c6_c8[1][i];
      Ry = c6_c8[2][i];
    }

    if (Ly && Ry) {
      num++;

      ddxyz(xyz[Ry], xyz[Ly], dd);
      vec_orth(dd, dz); /* corrected y-axis */
      cpxyz(dd, dy);
      cross(dy, dz, dx);

      ddxyz(o1, o2, dd);
      pars[1] = dot(dd, dx);
      pars[2] = dot(dd, dy);
      pars[3] = dot(dd, dz);

      buckle = vec_ang(z1, z2, dx);
      propeller = vec_ang(z1, z2, dy);

      norm = sqrt(buckle * buckle + propeller * propeller);
      if (norm > 1.0e-6) {
        r = angle / norm;
        pars[4] = r * buckle;
        pars[5] = r * propeller;
      } else { /* ~0 */
        pars[4] = buckle;
        pars[5] = propeller;
      }

      pars[6] = vec_ang(y1, y2, dz);

      for (j = 1; j <= 6; j++) {
        sprintf(tmp, fmt, pars[j]);
        strcat(str, tmp);
        parcln[num][j] = pars[j];
      }

    } else {
      for (j = 1; j <= 6; j++)
        strcat(str, bstr);
    }

    fprintf(fp, "%s\n", str);
  }

  output_ave_std(num, parcln, 1, fmt, fp);

  free_dmatrix(r1, 1, DUMMY, 1, DUMMY);
  free_dmatrix(r2, 1, DUMMY, 1, DUMMY);
  free_dmatrix(mst, 1, DUMMY, 1, DUMMY);
  free_dmatrix(parcln, 1, DUMMY, 1, DUMMY);
}

static long C1C1_based_frame(long idx, long **chi, double **xyz,
                             double **rotmat) {
  double x_axis[4], y_axis[4], z_axis[4];
  long k, c1a, c1b;

  k = (idx - 1) * 4;
  c1a = chi[1][k + 2];
  c1b = chi[2][k + 2];

  if (!c1a || !c1b)
    return false;

  /* with the C1'--C1' vector, for y-axis */
  ddxyz(xyz[c1b], xyz[c1a], y_axis); /* 2-->1 */

  /* keep origin and z-axis unchanged */
  for (k = 1; k <= 3; k++)
    z_axis[k] = rotmat[k][3];

  vec_orth(y_axis, z_axis); /* corrected y-axis, orthogonal to the z-axis */
  cross(y_axis, z_axis, x_axis); /* x-axis */

  for (k = 1; k <= 3; k++) {
    rotmat[k][1] = x_axis[k];
    rotmat[k][2] = y_axis[k];
  }

  return true;
}

static void simple_step_heli_pars(double **orien, double **org, long **chi,
                                  long num_bp, char **bp_seq, long *WC_info,
                                  long *bphlx, double **xyz, long hel_pars,
                                  FILE *fp) {
  char wc, str[BUF512], tmp[BUF32];
  char *bstr = "      ----", *fmt = "%10.2f";
  long i, j, ip1, num = 0;
  double morg[4], o1[4], o2[4], dd[4], pars[7];
  double *p1, *p2, **r1, **r2, **mst, **Rotmat;
  double **parcln;

  print_sep(fp, '-', 76);
  if (hel_pars) {
    fprintf(fp, "Simple base-pair helical parameters based on consecutive "
                "C1'-C1' vectors\n");
    fprintf(fp, "      step       X-disp    Y-disp   h-Rise     Incl.       "
                "Tip   h-Twist\n");
  } else {
    fprintf(fp, "Simple base-pair step parameters based on consecutive C1'-C1' "
                "vectors\n");
    fprintf(fp, "      step       Shift     Slide      Rise      Tilt      "
                "Roll     Twist\n");
  }

  r1 = dmatrix(1, 3, 1, 3);
  r2 = dmatrix(1, 3, 1, 3);
  mst = dmatrix(1, 3, 1, 3);
  Rotmat = dmatrix(1, num_bp, 1, 12);

  /* collect original pair frames */
  for (i = 1; i <= num_bp; i++) {
    refs_right_left(i, orien, org, r1, o1, r2, o2);
    bpstep_par(r1, o1, r2, o2, pars, mst, morg);

    p1 = Rotmat[i];        /* as a shorthand */
    mst2orien(p1, 0, mst); /* keep z-axis as is */
    cpxyz(morg, p1 + 9);   /* keep bp origin in columns 10-12 */
  }

  parcln = dmatrix(1, num_bp, 1, 6);
  for (i = 1; i < num_bp; i++) {
    ip1 = i + 1;

    wc = (WC_info[i] == 2 && WC_info[ip1] == 2) ? ' ' : '*';
    sprintf(str, "%c %4ld %c%c/%c%c", wc, i, bp_seq[1][i], bp_seq[1][ip1],
            bp_seq[2][ip1], bp_seq[2][i]);

    p1 = Rotmat[i];
    p2 = Rotmat[ip1];
    orien2mst(p1, 0, r1);
    orien2mst(p2, 0, r2);

    /* make frames with their z-axes along org1--->org2 */
    ddxyz(p1 + 9, p2 + 9, dd);
    if (dot(dd, p1 + 6) < 0)
      reverse_y_z_columns(r1);
    if (dot(dd, p2 + 6) < 0)
      reverse_y_z_columns(r2);

    if (C1C1_based_frame(i, chi, xyz, r1) &&
        C1C1_based_frame(ip1, chi, xyz, r2) && !bphlx[i]) {
      num++;
      if (hel_pars)
        helical_par(r1, p1 + 9, r2, p2 + 9, pars, mst, morg);
      else
        bpstep_par(r1, p1 + 9, r2, p2 + 9, pars, mst, morg);
      for (j = 1; j <= 6; j++) {
        sprintf(tmp, fmt, pars[j]);
        strcat(str, tmp);
        parcln[num][j] = pars[j];
      }

    } else
      for (j = 1; j <= 6; j++)
        strcat(str, bstr);

    fprintf(fp, "%s\n", str);
  }

  output_ave_std(num, parcln, 2, fmt, fp);

  free_dmatrix(r1, 1, DUMMY, 1, DUMMY);
  free_dmatrix(r2, 1, DUMMY, 1, DUMMY);
  free_dmatrix(mst, 1, DUMMY, 1, DUMMY);
  free_dmatrix(Rotmat, 1, DUMMY, 1, DUMMY);
  free_dmatrix(parcln, 1, DUMMY, 1, DUMMY);
}

static void check_simple_parameters(long simple, double **orien, double **org,
                                    long **c6_c8, long **chi, long num_bp,
                                    char **bp_seq, long *WC_info, long *bphlx,
                                    double **xyz, FILE *fp) {
  long i, num, NNy, k = 0;

  if (!simple)
    return;

  for (i = 1; i <= num_bp; i++)
    if (WC_info[i] == 2)
      k++;
  num = num_bp - k; /* number of non-WC base-pairs */

  print_sep(fp, '*', 76);
  fprintf(fp, "The 'simple' parameters are intuitive for non-Watson-Crick base "
              "pairs\n");
  fprintf(
      fp,
      "and associated base-pair steps (where the above corresponding 3DNA\n");
  fprintf(
      fp,
      "parameters often appear cryptic). Note that they are for structural\n");
  fprintf(fp, "*description* only, not to be fed into the 'rebuild' program. "
              "See URL\n");
  fprintf(fp, "http://x3dna.org/highlights/"
              "details-on-the-simple-base-pair-parameters\n");
  fprintf(fp, "and related blogposts on the 3DNA home page for details.\n\n");

  fprintf(fp,
          "This structure contains %ld non-Watson-Crick (with leading *) base "
          "pair(s)\n",
          num);

  NNy = simple & SIMPLE_BP_LONG_AXIS_RN9_YN1; /* RN9--YN1 as (long) y-axis */
  simple_bp_pars(orien, org, NNy, c6_c8, chi, num_bp, bp_seq, WC_info, xyz, fp);

  simple_step_heli_pars(orien, org, chi, num_bp, bp_seq, WC_info, bphlx, xyz,
                        false, fp);
  if (simple &
      SIMPLE_STEP_HELICAL_PARS) /* also derive/output helical parameters */
    simple_step_heli_pars(orien, org, chi, num_bp, bp_seq, WC_info, bphlx, xyz,
                          true, fp);
}

static void process_str(char *inpfile, struct_args *args) {
  char pdbfile[BUF512], outfile[BUF512], **nt_info;
  char *ChainID, *bseq, **AtomName, **bp_seq, **ResName, **Miscs;
  double twist_p, twist_n, *mst_org, *mst_orgH, *mst_orien, *mst_orienH;
  double **org, **orien, **twist_rise, **xyz, **nt_bb_torsion;
  long i, ip, bbexist = 0, ds, hetatm = 0, parallel = 0, nbpm1, num;
  long num_bp, num_residue, str_type = 0; /* for duplex only */
  long *idx, *WC_info, *ResSeq, *bphlx, *RY;
  long **c6_c8, **chi, **pair_num, **phos;
  long **seidx, **sugar, **o3p_brk, **htm_water;
  FILE *fp;

  /* get PDB file and pairing information */
  pair_num = read_input(inpfile, pdbfile, outfile, &ds, &num_bp, &ip, &hetatm);

  bphlx = lvector(1, num_bp); /* helix break marker */
  for (i = 1; i <= num_bp; i++)
    bphlx[i] = args->icnt ? 0 : pair_num[ds + 1][i];

  fp = open_file(outfile, "w");

  /* read in PDB file */
  num = number_of_atoms(pdbfile, hetatm, Gvars.misc_pars.alt_list);
  AtomName = cmatrix(1, num, 0, 4);
  ResName = cmatrix(1, num, 0, 3);
  ChainID = cvector(1, num);
  ResSeq = lvector(1, num);
  xyz = dmatrix(1, num, 1, 3);
  Miscs = cmatrix(1, num, 0, NMISC);
  read_pdb(pdbfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
           hetatm, Gvars.misc_pars.alt_list);

  idx = lvector(1, num);
  atom_idx(num, AtomName, NULL, idx);

  /* get the numbering information of each residue */
  seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

  /* base-pairing residue number checking */
  pair_checking(ip, ds, num_residue, pdbfile, &num_bp, pair_num);

  /* check strand direction: O3'-P linkage, parallel etc */
  o3p_brk = lmatrix(1, ds, 1, num_bp);
  drct_checking(ds, num_bp, pair_num, seidx, AtomName, xyz, &parallel, &bbexist,
                o3p_brk, fp);

  /* get base or base-pair sequence */
  bp_seq = cmatrix(0, ds, 1, num_bp); /* need to be kept */
  RY = lvector(1, num_residue);       /* simplified */
  get_bpseq(ds, num_bp, pair_num, seidx, AtomName, ResName, ChainID, ResSeq,
            Miscs, xyz, bp_seq, RY);

  /* for base overlap */
  bseq = cvector(1, num_residue);
  get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
          bseq, RY); /* RY re-calculated */

  print_header(ds, num_bp, num, pdbfile, fp);

  /* get atom list for each residue */
  phos = lmatrix(1, ds + 4, 1, num_bp); /* also include O1P & O2P atoms */
  c6_c8 = lmatrix(1, ds, 1, num_bp);
  sugar = lmatrix(1, ds, 1, num_bp * 5);
  chi = lmatrix(1, ds, 1, num_bp * 4);
  atom_list(ds, num_bp, pair_num, seidx, RY, bp_seq, AtomName, ResName, ChainID,
            ResSeq, Miscs, phos, c6_c8, sugar, chi);

  nt_info = cmatrix(1, num_residue, 0, BUF32);
  populate_nt_info(num_residue, seidx, ResName, ChainID, ResSeq, Miscs, bseq,
                   nt_info);

  /* get the reference frame for each base */
  org = dmatrix(1, ds, 1, num_bp * 3);
  orien = dmatrix(1, ds, 1, num_bp * 9);
  WC_info = lvector(1, num_bp);
  ref_frames(ds, num_bp, pair_num, bp_seq, seidx, RY, AtomName, ResName,
             ChainID, ResSeq, Miscs, xyz, fp, orien, org, WC_info, &str_type, 0,
             o3p_brk);

  nt_bb_torsion = dmatrix(1, num_residue, 1, 6);
  get_nt_bb_torsion(nt_bb_torsion, num_residue, seidx, RY, AtomName, ResName,
                    ChainID, ResSeq, Miscs, xyz);

  if (ds == 2) /* H-Bond information */
    hb_information(num_bp, pair_num, bp_seq, seidx, idx, AtomName, xyz, WC_info,
                   fp);

  base_overlap(ds, num_bp, num, num_residue, pair_num, RY, bp_seq, seidx,
               AtomName, xyz, idx, orien, org,
               fp); /* base overlap area in A^2 */

  /* get and print out parameters */
  nbpm1 = num_bp - 1;
  mst_orien = dvector(1, nbpm1 * 9);
  mst_org = dvector(1, nbpm1 * 3);
  mst_orienH = dvector(1, nbpm1 * 9);
  mst_orgH = dvector(1, nbpm1 * 3);
  twist_rise = dmatrix(1, nbpm1, 1, 2);
  get_parameters(ds, num_bp, bp_seq, orien, org, WC_info, fp, twist_rise,
                 mst_orien, mst_org, mst_orienH, mst_orgH, bphlx, args->istart,
                 args->istep, args->bz, &str_type, pair_num, nt_info);

  htm_water = lmatrix(1, 4, 0, num); /* HETATM and water index */
  init_htm_water(args->waters, num, num_residue, idx, htm_water);
  identify_htw(num_residue, seidx, RY, AtomName, ResName, ChainID, ResSeq,
               Miscs, xyz, htm_water);

  write_mst(ds, num_bp, pair_num, bp_seq, mst_orien, mst_org, seidx, AtomName,
            ResName, ChainID, ResSeq, xyz, Miscs, htm_water, twist_rise,
            STACK_FILE);
  write_mst(ds, num_bp, pair_num, bp_seq, mst_orienH, mst_orgH, seidx, AtomName,
            ResName, ChainID, ResSeq, xyz, Miscs, htm_water, twist_rise,
            HSTACK_FILE);

  if (args->ring)
    get_ring_center_normal(ds, num_bp, pair_num, bp_seq, seidx, RY, AtomName,
                           ResName, ChainID, ResSeq, Miscs, xyz, fp);

  if (ds == 2) {
    check_simple_parameters(args->simple_pars, orien, org, c6_c8, chi, num_bp,
                            bp_seq, WC_info, bphlx, xyz, fp);

    get_mtwist(nbpm1, bphlx, WC_info, twist_rise, &twist_p, &twist_n);

    str_classify(twist_p, twist_n, str_type, parallel, num_bp, fp);
    lambda_d3(num_bp, bp_seq, chi, c6_c8, xyz, fp);

    if (bbexist) {
      print_PP(parallel, twist_rise, num_bp, bp_seq, phos, mst_orien, mst_org,
               mst_orienH, mst_orgH, xyz, WC_info, bphlx, args->abi, chi, fp);
      groove_width(parallel, num_bp, bp_seq, phos, xyz, bphlx, fp);
    }
    other_pars(num_bp, bp_seq, bphlx, orien, org);
  }

  /* global analysis */
  global_analysis(ds, num_bp, num, bp_seq, chi, phos, xyz, fp);

  if (bbexist) {
    backbone_torsion(ds, num_bp, pair_num, bp_seq, sugar, chi, xyz,
                     nt_bb_torsion, fp);
    p_c1_dist(ds, num_bp, bp_seq, phos, chi, xyz, bphlx, fp);
    helix_radius(ds, num_bp, bp_seq, orien, org, phos, chi, xyz, bphlx, fp);
  }
  get_helix_axis(ds, num_bp, bp_seq, orien, org, bphlx, fp);

  /* free allocated vectors & matrices */
  free_lmatrix(pair_num, 1, ds + 1, 1, num_bp);
  free_lvector(bphlx, 1, num_bp);
  free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
  free_lmatrix(seidx, 1, num_residue, 1, 2);
  free_lmatrix(o3p_brk, 1, ds, 1, num_bp);
  free_cmatrix(bp_seq, 0, ds, 1, num_bp);
  free_lvector(RY, 1, num_residue);
  free_cvector(bseq, 1, num_residue);
  free_lmatrix(phos, 1, ds + 4, 1, num_bp);
  free_lmatrix(c6_c8, 1, ds, 1, num_bp);
  free_dmatrix(nt_bb_torsion, 1, num_residue, 1, 6);
  free_lmatrix(sugar, 1, ds, 1, num_bp * 5);
  free_lmatrix(chi, 1, ds, 1, num_bp * 4);
  free_dmatrix(orien, 1, ds, 1, num_bp * 9);
  free_dmatrix(org, 1, ds, 1, num_bp * 3);
  free_lvector(WC_info, 1, num_bp);
  free_dvector(mst_orien, 1, nbpm1 * 9);
  free_dvector(mst_org, 1, nbpm1 * 3);
  free_dvector(mst_orienH, 1, nbpm1 * 9);
  free_dvector(mst_orgH, 1, nbpm1 * 3);
  free_dmatrix(twist_rise, 1, nbpm1, 1, 2);
  free_lmatrix(htm_water, 1, 4, 0, num);
  free_lvector(idx, 1, num);
  free_cmatrix(nt_info, 1, num_residue, 0, BUF32);

  close_file(fp);
}

static void calculate_torsions(char *outfile, char *pdbfile) {
  char BDIR[BUF512], **nt_info;
  char *ChainID, *bseq, **AtomName, **ResName, **Miscs;
  double **org, **orien, **xyz, **nt_torsion, **ss_Zp_Dp;
  long num, num_residue, hetatm = true;
  long *ResSeq, *RY, **seidx, **nt_list;
  FILE *fp;

  fp = open_file(outfile, "w");

  num = number_of_atoms(pdbfile, hetatm, Gvars.misc_pars.alt_list);
  AtomName = cmatrix(1, num, 0, 4);
  ResName = cmatrix(1, num, 0, 3);
  ChainID = cvector(1, num);
  ResSeq = lvector(1, num);
  xyz = dmatrix(1, num, 1, 3);
  Miscs = cmatrix(1, num, 0, NMISC);
  read_pdb(pdbfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
           hetatm, Gvars.misc_pars.alt_list);

  seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

  bseq = cvector(1, num_residue);
  RY = lvector(1, num_residue);
  get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
          bseq, RY);

  org = dmatrix(1, num_residue, 1, 3);
  orien = dmatrix(1, num_residue, 1, 9);
  get_BDIR(BDIR, "Atomic_A.pdb");
  base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq,
             Miscs, xyz, BDIR, orien, org);

  nt_info = cmatrix(1, num_residue, 0, BUF32);
  populate_nt_info(num_residue, seidx, ResName, ChainID, ResSeq, Miscs, bseq,
                   nt_info);

  nt_list = lmatrix(1, num_residue, 0, BUF32);
  populate_nt_list(num_residue, seidx, RY, bseq, AtomName, xyz, nt_list);
  output_Borg_P_C1_C4(num_residue, org, xyz, nt_list, nt_info);

  nt_torsion = dmatrix(1, num_residue, 1, BUF32);
  get_nt_torsion(num_residue, org, xyz, nt_list, nt_torsion);

  ss_Zp_Dp = dmatrix(1, 2, 1, num_residue);
  get_ss_Zp_Dp(num_residue, org, orien, xyz, nt_list, ss_Zp_Dp);

  output_nt_torsion(num_residue, nt_info, nt_list, nt_torsion, ss_Zp_Dp, fp);

  free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
  free_lmatrix(seidx, 1, num_residue, 1, 2);
  free_lvector(RY, 1, num_residue);
  free_cvector(bseq, 1, num_residue);
  free_dmatrix(org, 1, num_residue, 1, 3);
  free_dmatrix(orien, 1, num_residue, 1, 9);
  free_dmatrix(nt_torsion, 1, num_residue, 1, BUF32);
  free_lmatrix(nt_list, 1, num_residue, 0, BUF32);
  free_cmatrix(nt_info, 1, num_residue, 0, BUF32);
  free_dmatrix(ss_Zp_Dp, 1, 2, 1, num_residue);

  close_file(fp);
}

static void analyze_usage(void) { help3dna_usage("analyze"); }

static void set_defaults(struct_args *args) {
  strcpy(args->torsion, "");
  args->istart = 1;
  args->istep = 1;
  args->icnt = false;
  args->waters = false;
  args->bz = true; /* check for B-Z junction */
  args->ring = false;
  args->simple_pars = true;
  args->abi = false;
  args->circular = false;
}

static long derive_simple_pars_settings(char *option) {
  long k = false;

  if (lux_ncmatch(option, "no|false|off"))
    return k;

  k = true; /* default: RC8--YC6, no helical pars */

  if (lux_ncmatch(option, "n1|n9"))
    k |= SIMPLE_BP_LONG_AXIS_RN9_YN1;

  if (lux_ncmatch(option, "heli"))
    k |= SIMPLE_STEP_HELICAL_PARS;

  fprintf(stderr, "Setting for simple parameters: %ld\n", k);
  return k;
}

static void analyze_cmdline(int argc, char *argv[], struct_args *args) {
  long i, j, k;

  set_defaults(args);

  for (i = 1; i < argc; i++) {
    if (*argv[i] != '-')
      break;

    if (check_global_options(argv[i]))
      continue;

    if (lux_ncmatch(argv[i], "^--?t")) {
      get_strvalue(argv[i], args->torsion, false);
      continue;
    }

    if (lux_ncmatch(argv[i], "^--?bz")) {
      args->bz = set_switch_default_true(argv[i]);
      continue;
    }

    if (lux_ncmatch(argv[i], "^--?ri")) {
      args->ring = set_switch_default_true(argv[i]);
      continue;
    }

    if (lux_ncmatch(argv[i], "^--?si")) {
      args->simple_pars = derive_simple_pars_settings(argv[i]);
      continue;
    }

    if (lux_ncmatch(argv[i], "^--?abi?")) {
      args->abi = set_switch_default_true(argv[i]);
      continue;
    }

    if (lux_ncmatch(argv[i], "^--?circ")) {
      args->circular = set_switch_default_true(argv[i]);
      continue;
    }

    upperstr(argv[i]);
    for (j = 1; j < (long)strlen(argv[i]); j++) /* skip - */
      if (argv[i][j] == 'C')
        args->icnt = true;
      else if (argv[i][j] == 'W')
        args->waters = true;
      else if (argv[i][j] == 'S') {
        k = sscanf(argv[i], "-S=%ld,%ld", &args->istep, &args->istart);
        if (k == 2) {
          if (args->istart < 0)
            args->istart = -args->istart; /* make it positive */
        } else if (k == 1)
          args->istart = 1; /* redundant: already initialized */
        else {
          fprintf(stderr, "wrong format for setting step\n");
          analyze_usage();
        }
        fprintf(stderr, "***start at %ld, with step size: %ld***\n",
                args->istart, args->istep);
        break; /* -s=istep,istart not combined with others */
      } else
        analyze_usage();
  }

  if (argc == i) {
    if (is_empty_string(args->torsion))
      process_str("stdin", args);
    else
      analyze_usage();
  } else {
    if (is_empty_string(args->torsion)) { /* regular case */
      for (j = i; j < argc; j++) {
        if (strcmp(argv[j], "tmpfile")) /* analyze multiple structures */
          fprintf(stderr, "\n......Processing structure #%ld: <%s>......\n",
                  j - i + 1, argv[j]);
        process_str(argv[j], args);
      }
    } else {                                      /* torsion angle only */
      calculate_torsions(args->torsion, argv[i]); /* only one PDB structure */
    }
  }
}
/*
int main(int argc, char *argv[]) {
  struct_args args;
  time_t time0;

  time(&time0);

  // clean up files with fixed name
  remove_file(AUX_FILE);
  remove_file(BPSTEP_FILE);
  remove_file(HLXSTEP_FILE);
  remove_file(STACK_FILE);
  remove_file(HSTACK_FILE);
  remove_file(REF_FILE);
  remove_file(POC_FILE);
  remove_file(SEVEN_FILE);

  set_my_globals(argv[0]);
  analyze_cmdline(argc, argv, &args);
  clear_my_globals();

  print_used_time(time0);

  return 0;
}
*/

/// @brief - nrutil.cpp

/* standard error handler */
void nrerror(char *error_text) {
  fprintf(stderr, "%s\n", error_text);
  exit(1);
}

void vector_boundary_check(long nl, long nh, char *fun_name) {
  return; /* not used */

  if (nl > nh)
    fprintf(stderr, "boundary for %s: low = %ld > high = %ld\n", fun_name, nl,
            nh);
}

void matrix_boundary_check(long nrl, long nrh, long ncl, long nch,
                           char *fun_name) {
  return; /* not used */

  if (nrl > nrh || ncl > nch)
    fprintf(stderr, "boundary for %s: [%ld to %ld; %ld to %ld]\n", fun_name,
            nrl, nrh, ncl, nch);
}

/* ------------------------------------------------------------------ */
/* allocate a char vector with subscript range v[nl..nh] */
char *cvector(long nl, long nh) {
  char *v;

  vector_boundary_check(nl, nh, "cvector()");
  v = cvector_nr(nl, nh);
  init_cvector(v, nl, nh, '\0');

  return v;
}

/* allocate a char vector with subscript range v[nl..nh] */
char *cvector_nr(long nl, long nh) {
  char *v;

  v = (char *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(char)));
  if (!v)
    nrerror("allocation failure in cvector()");

  return v - nl + NR_END;
}

/* ------------------------------------------------------------------ */
/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(long nl, long nh) {
  double *v;

  vector_boundary_check(nl, nh, "dvector()");
  v = dvector_nr(nl, nh);
  init_dvector(v, nl, nh, 0.0);

  return v;
}

/* allocate a double vector with subscript range v[nl..nh] */
double *dvector_nr(long nl, long nh) {
  double *v;

  v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
  if (!v)
    nrerror("allocation failure in dvector()");

  return v - nl + NR_END;
}

/* ------------------------------------------------------------------ */
/* allocate a long vector with subscript range v[nl..nh] */
long *lvector(long nl, long nh) {
  long *v;

  vector_boundary_check(nl, nh, "lvector()");
  v = lvector_nr(nl, nh);
  init_lvector(v, nl, nh, 0);

  return v;
}

/* allocate a long vector with subscript range v[nl..nh] */
long *lvector_nr(long nl, long nh) {
  long *v;

  v = (long *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(long)));
  if (!v)
    nrerror("allocation failure in lvector()");

  return v - nl + NR_END;
}

/* ------------------------------------------------------------------ */
/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */

char **cmatrix(long nrl, long nrh, long ncl, long nch) {
  char **m;

  matrix_boundary_check(nrl, nrh, ncl, nch, "cmatrix()");
  m = cmatrix_nr(nrl, nrh, ncl, nch);
  init_cmatrix(m, nrl, nrh, ncl, nch, '\0');

  return m;
}

/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
char **cmatrix_nr(long nrl, long nrh, long ncl, long nch) {
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  char **m;

  /* allocate pointers to rows */
  m = (char **)malloc((size_t)((nrow + NR_END) * sizeof(char *)));
  if (!m)
    nrerror("allocation failure 1 in cmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (char *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(char)));
  if (!m[nrl])
    nrerror("allocation failure 2 in cmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* ------------------------------------------------------------------ */
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(long nrl, long nrh, long ncl, long nch) {
  double **m;

  matrix_boundary_check(nrl, nrh, ncl, nch, "dmatrix()");
  m = dmatrix_nr(nrl, nrh, ncl, nch);
  init_dmatrix(m, nrl, nrh, ncl, nch, 0.0);

  return m;
}

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix_nr(long nrl, long nrh, long ncl, long nch) {
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double *)));
  if (!m)
    nrerror("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
  if (!m[nrl])
    nrerror("allocation failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* ------------------------------------------------------------------ */
/* allocate a long matrix with subscript range m[nrl..nrh][ncl..nch] */
long **lmatrix(long nrl, long nrh, long ncl, long nch) {
  long **m;

  matrix_boundary_check(nrl, nrh, ncl, nch, "lmatrix()");
  m = lmatrix_nr(nrl, nrh, ncl, nch);
  init_lmatrix(m, nrl, nrh, ncl, nch, 0);

  return m;
}

/* allocate a long matrix with subscript range m[nrl..nrh][ncl..nch] */
long **lmatrix_nr(long nrl, long nrh, long ncl, long nch) {
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  long **m;

  /* allocate pointers to rows */
  m = (long **)malloc((size_t)((nrow + NR_END) * sizeof(long *)));
  if (!m)
    nrerror("allocation failure 1 in lmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (long *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(long)));
  if (!m[nrl])
    nrerror("allocation failure 2 in lmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* ------------------------------------------------------------------ */
/* free a char vector allocated with cvector() */
void free_cvector(char *v, long nl, long nh) {
  free((FREE_ARG)(v + nl - NR_END));
  UNUSED_PARAMETER(nh); /* strict compiler options */
}

/* free a double vector allocated with dvector() */
void free_dvector(double *v, long nl, long nh) {
  free((FREE_ARG)(v + nl - NR_END));
  UNUSED_PARAMETER(nh); /* strict compiler options */
}

/* free a long vector allocated with lvector() */
void free_lvector(long *v, long nl, long nh) {
  free((FREE_ARG)(v + nl - NR_END));
  UNUSED_PARAMETER(nh); /* strict compiler options */
}

/* free a char matrix allocated by cmatrix() */
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
  UNUSED_PARAMETER(nrh);
  UNUSED_PARAMETER(nch); /* strict compiler options */
}

/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
  UNUSED_PARAMETER(nrh);
  UNUSED_PARAMETER(nch); /* strict compiler options */
}

/* free a long matrix allocated by lmatrix() */
void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
  UNUSED_PARAMETER(nrh);
  UNUSED_PARAMETER(nch); /* strict compiler options */
}

/* ------------------------------------------------------------------ */
void dval_swap(double *pa, double *pb) {
  double temp;

  temp = *pa;
  *pa = *pb;
  *pb = temp;
}

void lval_swap(long *pa, long *pb) {
  long temp;

  temp = *pa;
  *pa = *pb;
  *pb = temp;
}

void cval_swap(char *pa, char *pb) {
  int c;

  c = *pa;
  *pa = *pb;
  *pb = c;
}

double dval_max(double a, double b) { return (a > b) ? a : b; }

double dval_min(double a, double b) { return (a < b) ? a : b; }

long lval_min(long a, long b) { return (a < b) ? a : b; }

/* check if lval exists in s[ib] .. s[ie] */
long lval_in_set(long lval, long ib, long ie, long *s) {
  long i;

  for (i = ib; i <= ie; i++)
    if (lval == s[i])
      return true;

  return false;
}

/* check if "dval" in within "dlow" and "dhigh" */
long dval_in_range(double dval, double dlow, double dhigh) {
  return (dval >= dlow && dval <= dhigh) ? 1L : 0L;
}

/* check if "lval" in within "llow" and "lhigh" */
long lval_in_range(long lval, long llow, long lhigh) {
  return (lval >= llow && lval <= lhigh) ? 1L : 0L;
}

void max_dmatrix(double **d, long nr, long nc, double *maxdm) {
  long i, j;

  for (i = 1; i <= nc; i++) {
    maxdm[i] = -XBIG;
    for (j = 1; j <= nr; j++)
      maxdm[i] = dval_max(maxdm[i], d[j][i]);
  }
}

void min_dmatrix(double **d, long nr, long nc, double *mindm) {
  long i, j;

  for (i = 1; i <= nc; i++) {
    mindm[i] = XBIG;
    for (j = 1; j <= nr; j++)
      mindm[i] = dval_min(mindm[i], d[j][i]);
  }
}

void ave_dmatrix(double **d, long nr, long nc, double *avedm) {
  long i, j;

  if (!nr)
    nrerror("divided by zero in <ave_dmatrix>");

  for (i = 1; i <= nc; i++) {
    avedm[i] = 0.0;
    for (j = 1; j <= nr; j++)
      avedm[i] += d[j][i];
    avedm[i] /= nr;
  }
}

void std_dmatrix(double **d, long nr, long nc, double *stddm) {
  double dsum, temp;
  double *aved;
  long i, j;

  if (nr < 2) {
    fprintf(stderr, "less than two items for <std_dmatrix>\n");
    init_dvector(stddm, 1, nc, 0.0); /* column-wise */
    return;
  }

  aved = dvector(1, nc);
  ave_dmatrix(d, nr, nc, aved);

  for (i = 1; i <= nc; i++) {
    dsum = 0.0;
    for (j = 1; j <= nr; j++) {
      temp = d[j][i] - aved[i];
      dsum += temp * temp;
    }
    stddm[i] = sqrt(dsum / (nr - 1));
  }

  free_dvector(aved, 1, nc);
}

double max_dvector(double *d, long nl, long nh) {
  double maxdv = -XBIG;
  long i;

  for (i = nl; i <= nh; i++)
    maxdv = dval_max(maxdv, d[i]);

  return maxdv;
}

double min_dvector(double *d, long nl, long nh) {
  double mindv = XBIG;
  long i;

  for (i = nl; i <= nh; i++)
    mindv = dval_min(mindv, d[i]);

  return mindv;
}

double ave_dvector(double *d, long n) {
  double dsum = 0.0;
  long i;

  if (n <= 0) {
    fprintf(stderr, "zero (or negative) items [%ld] for <ave_dvector>\n", n);
    return 0.0;
  }

  for (i = 1; i <= n; i++)
    dsum += d[i];

  return dsum / n;
}

double std_dvector(double *d, long n) {
  double aved, dsum = 0.0, temp;
  long i;

  if (n < 2) {
    fprintf(stderr, "less than 2 items for <std_dvector>\n");
    return 0.0;
  }

  aved = ave_dvector(d, n);
  for (i = 1; i <= n; i++) {
    temp = d[i] - aved;
    dsum += temp * temp;
  }

  return sqrt(dsum / (n - 1));
}

/* initialize character-matrix dmtx with init_val */
void init_cmatrix(char **cmtx, long nrl, long nrh, long ncl, long nch,
                  char init_val) {
  long i;

  for (i = nrl; i <= nrh; i++)
    init_cvector(cmtx[i], ncl, nch, init_val);
}

/* initialize double-matrix dmtx with init_val */
void init_dmatrix(double **dmtx, long nrl, long nrh, long ncl, long nch,
                  double init_val) {
  long i;

  for (i = nrl; i <= nrh; i++)
    init_dvector(dmtx[i], ncl, nch, init_val);
}

/* initialize long-matrix lmtx with init_val */
void init_lmatrix(long **lmtx, long nrl, long nrh, long ncl, long nch,
                  long init_val) {
  long i;

  for (i = nrl; i <= nrh; i++)
    init_lvector(lmtx[i], ncl, nch, init_val);
}

/* initialize a character vector cvec[ib..ie] with init_val, as a string */
void init_cvector(char *cvec, long ib, long ie, char init_val) {
  long i;

  for (i = ib; i < ie; i++)
    cvec[i] = init_val;
  cvec[ie] = '\0';
}

/* initialize a double vector dvec[ib..ie] with init_val */
void init_dvector(double *dvec, long ib, long ie, double init_val) {
  long i;

  for (i = ib; i <= ie; i++)
    dvec[i] = init_val;
}

/* initialize a long vector lvec[ib..ie] with init_val */
void init_lvector(long *lvec, long ib, long ie, long init_val) {
  long i;

  for (i = ib; i <= ie; i++)
    lvec[i] = init_val;
}

/* copy a double vector [nl .. nh], from s to d */
void copy_dvector(double *d, double *s, long nl, long nh) {
  long i;

  for (i = nl; i <= nh; i++)
    d[i] = s[i];
}

int dval_compare(const void *v1, const void *v2) {
  const double *p1, *p2;

  p1 = (const double *)v1;
  p2 = (const double *)v2;

  return (*p1 > *p2) ? 1 : (*p1 < *p2) ? -1 : 0;
}

/* negate each value of vector xyz1 */
void negate_xyz(double *xyz1) {
  long i;

  for (i = 1; i <= 3; i++)
    xyz1[i] = -xyz1[i];
}

/* distance between vectors xyz1 and xyz2 */
double p1p2_dist(double *xyz1, double *xyz2) {
  double dxyz[4];

  ddxyz(xyz1, xyz2, dxyz);
  return veclen(dxyz);
}

/* check if distance between xyz1 & xyz2 is within (dlow, dhigh) */
long within_limits(double *xyz1, double *xyz2, double dlow, double dhigh) {
  double d, dxyz[4];
  long i, ltok = 0;

  for (i = 1; i <= 3; i++)
    if (fabs(dxyz[i] = xyz1[i] - xyz2[i]) > dhigh)
      break;
  if (i > 3) {
    d = veclen(dxyz);
    if (dval_in_range(d, dlow, dhigh))
      ltok = 1;
  }
  return ltok;
}

/* get the sum vector of xyz1 and xyz2 */
void sumxyz(double *xyz1, double *xyz2, double *sxyz) {
  long i;

  for (i = 1; i <= 3; i++)
    sxyz[i] = xyz1[i] + xyz2[i];
}

/* get the mean vector of xyz1 and xyz2 */
void avexyz(double *xyz1, double *xyz2, double *mxyz) {
  long i;

  for (i = 1; i <= 3; i++)
    mxyz[i] = 0.5 * (xyz1[i] + xyz2[i]);
}

/* the difference vector of xyz2 and xyz1 */
void ddxyz(double *xyz1, double *xyz2, double *dxyz) {
  long i;

  for (i = 1; i <= 3; i++)
    dxyz[i] = xyz2[i] - xyz1[i];
}

/* make a copy of xyz1 to xyz2 */
void cpxyz(double *xyz1, double *xyz2) {
  long i;

  for (i = 1; i <= 3; i++)
    xyz2[i] = xyz1[i];
}

/* get orthogonal component of va w.r.t. vref [1-by-3] */
void vec_orth(double *va, double *vref) {
  double d;
  long i;

  vec_norm(vref);
  d = dot(va, vref);
  for (i = 1; i <= 3; i++)
    va[i] -= d * vref[i];
  vec_norm(va);
}

/* dot product between two 1-by-3 vectors */
double dot(double *va, double *vb) {
  double dsum = 0.0;
  long i;

  for (i = 1; i <= 3; i++)
    dsum += va[i] * vb[i];

  return dsum;
}

/* cross product between two 1-by-3 vectors */
void cross(double *va, double *vb, double *vc) {
  vc[1] = va[2] * vb[3] - va[3] * vb[2];
  vc[2] = va[3] * vb[1] - va[1] * vb[3];
  vc[3] = va[1] * vb[2] - va[2] * vb[1];
}

long sign_control(double *va, double *vb, double *vref) {
  double dtmp[4];

  cross(va, vb, dtmp);
  if (dot(dtmp, vref) < 0.0)
    return -1;
  else
    return 1;
}

/* length (magnitude) of a 1-by-3 vector */
double veclen(double *va) { return sqrt(dot(va, va)); }

/* normalize a 1-by-3 vector, i.e. with unit length */
void vec_norm(double *va) {
  double vlen;
  long i;

  vlen = veclen(va);
  if (vlen > XEPS)
    for (i = 1; i <= 3; i++)
      va[i] /= vlen;
}

/* from dot product of 2 NORMAL vectors to angle: positive */
double dot2ang(double dotval) {
  double ang_deg;
  if (dotval >= 1.0)
    ang_deg = 0.0;
  else if (dotval <= -1.0)
    ang_deg = 180.0;
  else
    ang_deg = rad2deg(acos(dotval));
  return ang_deg;
}

/* angle magnitude in degrees between two 1-by-3 vectors */
double magang(double *va, double *vb) {
  double ang_deg;
  if (veclen(va) < XEPS || veclen(vb) < XEPS)
    ang_deg = 0.0;
  else {
    vec_norm(va);
    vec_norm(vb);
    ang_deg = dot2ang(dot(va, vb));
  }
  return ang_deg;
}

double rad2deg(double ang) { return ang * 180.0 / PI; }

double deg2rad(double ang) { return ang * PI / 180.0; }

void copy_dmatrix(double **a, long nr, long nc, double **o) {
  long i, j;

  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
      o[i][j] = a[i][j];
}

void multi_matrix(double **a, long nra, long nca, double **b, long nrb,
                  long ncb, double **o) {
  long i, j, k;

  if (nca != nrb)
    nrerror("matrices a and b do not conform");

  for (i = 1; i <= nra; i++) {
    for (j = 1; j <= ncb; j++) {
      o[i][j] = 0.0;
      for (k = 1; k <= nca; k++)
        o[i][j] += a[i][k] * b[k][j];
    }
  }
}

/* vector-matrix multiplication */
void multi_vec_matrix(double *a, long n, double **b, long nr, long nc,
                      double *o) {
  long i, j;

  if (n != nr)
    nrerror("vector and matrix do not conform");

  for (i = 1; i <= nc; i++) {
    o[i] = 0.0;
    for (j = 1; j <= n; j++)
      o[i] += a[j] * b[j][i];
  }
}

/* vector - transpose-of-matrix multiplication */
void multi_vec_Tmatrix(double *a, long n, double **b, long nr, long nc,
                       double *o) {
  double **tb; /* transpose of b */

  tb = dmatrix(1, nc, 1, nr);

  transpose_matrix(b, nr, nc, tb);
  multi_vec_matrix(a, n, tb, nc, nr, o);

  free_dmatrix(tb, 1, nc, 1, nr);
}

void transpose_matrix(double **a, long nr, long nc, double **o) {
  long i, j;

  for (i = 1; i <= nc; i++)
    for (j = 1; j <= nr; j++)
      o[i][j] = a[j][i];
}

void identity_matrix(double **d, long n) {
  long i;

  for (i = 1; i <= n; i++) {
    init_dvector(d[i], 1, n, 0.0);
    d[i][i] = 1.0;
  }
}

/// @brief - fncs_slre.cpp

/*  Copyright (c) 2004-2012 Sergey Lyubka <valenok@gmail.com> */
/*  All rights reserved */

/*  Permission is hereby granted, free of charge, to any person obtaining a copy
 */
/*  of this software and associated documentation files (the "Software"), to
 * deal */
/*  in the Software without restriction, including without limitation the rights
 */
/*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
/*  copies of the Software, and to permit persons to whom the Software is */
/*  furnished to do so, subject to the following conditions: */

/*  The above copyright notice and this permission notice shall be included in
 */
/*  all copies or substantial portions of the Software. */

/*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 */
/*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
/*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 */
/*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
/*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, */
/*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN */
/*  THE SOFTWARE. */

struct slre {
  unsigned char code[256];
  unsigned char data[256];
  int code_size;
  int data_size;
  int num_caps;
  int anchored;
  enum slre_option options;
  const char *error_string;
};

struct cap {
  const char *ptr;
  int len;
};

enum {
  END,
  BRANCH,
  ANY,
  EXACT,
  ANYOF,
  ANYBUT,
  OPEN,
  CLOSE,
  BOL,
  EOL,
  STAR,
  PLUS,
  STARQ,
  PLUSQ,
  QUEST,
  SPACE,
  NONSPACE,
  DIGIT
};

static const char *meta_characters = "|.^$*+?()[\\";
static const char *error_no_match = "No match";

static void set_jump_offset(struct slre *r, int pc, int offset) {
  assert(offset < r->code_size);
  if (r->code_size - offset > 0xff) {
    r->error_string = "Jump offset is too big";
  } else {
    r->code[pc] = (unsigned char)(r->code_size - offset);
  }
}

static void emit(struct slre *r, int code) {
  if (r->code_size >= (int)NELEMS(r->code)) {
    r->error_string = "RE is too long (code overflow)";
  } else {
    r->code[r->code_size++] = (unsigned char)code;
  }
}

static void store_char_in_data(struct slre *r, int ch) {
  if (r->data_size >= (int)sizeof(r->data)) {
    r->error_string = "RE is too long (data overflow)";
  } else {
    r->data[r->data_size++] = ch;
  }
}

static void exact(struct slre *r, const char **re) {
  int old_data_size = r->data_size;

  while (**re != '\0' && (strchr(meta_characters, **re)) == NULL) {
    store_char_in_data(r, *(*re)++);
  }

  emit(r, EXACT);
  emit(r, old_data_size);
  emit(r, r->data_size - old_data_size);
}

static int get_escape_char(const char **re) {
  int res;

  switch (*(*re)++) {
  case 'n':
    res = '\n';
    break;
  case 'r':
    res = '\r';
    break;
  case 't':
    res = '\t';
    break;
  case '0':
    res = 0;
    break;
  case 'S':
    res = NONSPACE << 8;
    break;
  case 's':
    res = SPACE << 8;
    break;
  case 'd':
    res = DIGIT << 8;
    break;
  default:
    res = (*re)[-1];
    break;
  }

  return res;
}

static void anyof(struct slre *r, const char **re) {
  int esc, old_data_size = r->data_size, op = ANYOF;

  if (**re == '^') {
    op = ANYBUT;
    (*re)++;
  }

  while (**re != '\0')
    switch (*(*re)++) {
    case ']':
      emit(r, op);
      emit(r, old_data_size);
      emit(r, r->data_size - old_data_size);
      return;

      break;
    case '\\':
      esc = get_escape_char(re);
      if ((esc & 0xff) == 0) {
        store_char_in_data(r, 0);
        store_char_in_data(r, esc >> 8);
      } else {
        store_char_in_data(r, esc);
      }
      break;
    default:
      store_char_in_data(r, (*re)[-1]);
      break;
    }

  r->error_string = "No closing ']' bracket";
}

static void relocate(struct slre *r, int begin, int shift) {
  emit(r, END);
  memmove(r->code + begin + shift, r->code + begin, r->code_size - begin);
  r->code_size += shift;
}

static void quantifier(struct slre *r, int prev, int op) {
  if (r->code[prev] == EXACT && r->code[prev + 2] > 1) {
    r->code[prev + 2]--;
    emit(r, EXACT);
    emit(r, r->code[prev + 1] + r->code[prev + 2]);
    emit(r, 1);
    prev = r->code_size - 3;
  }
  relocate(r, prev, 2);
  r->code[prev] = op;
  set_jump_offset(r, prev + 1, prev);
}

static void exact_one_char(struct slre *r, int ch) {
  emit(r, EXACT);
  emit(r, r->data_size);
  emit(r, 1);
  store_char_in_data(r, ch);
}

static void fixup_branch(struct slre *r, int fixup) {
  if (fixup > 0) {
    emit(r, END);
    set_jump_offset(r, fixup, fixup - 2);
  }
}

static void compile(struct slre *r, const char **re) {
  int op, esc, branch_start, last_op, fixup, cap_no, level;

  fixup = 0;
  level = r->num_caps;
  branch_start = last_op = r->code_size;

  for (;;)
    switch (*(*re)++) {

    case '\0':
      (*re)--;
      return;

      break;

    case '^':
      emit(r, BOL);
      break;

    case '$':
      emit(r, EOL);
      break;

    case '.':
      last_op = r->code_size;
      emit(r, ANY);
      break;

    case '[':
      last_op = r->code_size;
      anyof(r, re);
      break;

    case '\\':
      last_op = r->code_size;
      esc = get_escape_char(re);
      if (esc & 0xff00) {
        emit(r, esc >> 8);
      } else {
        exact_one_char(r, esc);
      }
      break;

    case '(':
      last_op = r->code_size;
      cap_no = ++r->num_caps;
      emit(r, OPEN);
      emit(r, cap_no);

      compile(r, re);
      if (*(*re)++ != ')') {
        r->error_string = "No closing bracket";
        return;
      }

      emit(r, CLOSE);
      emit(r, cap_no);
      break;

    case ')':
      (*re)--;
      fixup_branch(r, fixup);
      if (level == 0) {
        r->error_string = "Unbalanced brackets";
        return;
      }
      return;

      break;

    case '+':
    case '*':
      op = (*re)[-1] == '*' ? STAR : PLUS;
      if (**re == '?') {
        (*re)++;
        op = op == STAR ? STARQ : PLUSQ;
      }
      quantifier(r, last_op, op);
      break;

    case '?':
      quantifier(r, last_op, QUEST);
      break;

    case '|':
      fixup_branch(r, fixup);
      relocate(r, branch_start, 3);
      r->code[branch_start] = BRANCH;
      set_jump_offset(r, branch_start + 1, branch_start);
      fixup = branch_start + 2;
      r->code[fixup] = 0xff;
      break;

    default:
      (*re)--;
      last_op = r->code_size;
      exact(r, re);
      break;
    }
}

static const char *compile2(struct slre *r, const char *re) {
  r->error_string = NULL;
  r->code_size = r->data_size = r->num_caps = r->anchored = 0;

  if (*re == '^') {
    r->anchored++;
  }

  emit(r, OPEN);
  emit(r, 0);

  while (*re != '\0') {
    compile(r, &re);
  }

  if (r->code[2] == BRANCH) {
    fixup_branch(r, 4);
  }

  emit(r, CLOSE);
  emit(r, 0);
  emit(r, END);

  return r->error_string;
}

static const char *match(const struct slre *, int, const char *, int, int *,
                         struct cap *);

static void loop_greedy(const struct slre *r, int pc, const char *s, int len,
                        int *ofs) {
  int saved_offset, matched_offset;

  saved_offset = matched_offset = *ofs;

  while (!match(r, pc + 2, s, len, ofs, NULL)) {
    saved_offset = *ofs;
    if (!match(r, pc + r->code[pc + 1], s, len, ofs, NULL)) {
      matched_offset = saved_offset;
    }
    *ofs = saved_offset;
  }

  *ofs = matched_offset;
}

static void loop_non_greedy(const struct slre *r, int pc, const char *s,
                            int len, int *ofs) {
  int saved_offset = *ofs;

  while (!match(r, pc + 2, s, len, ofs, NULL)) {
    saved_offset = *ofs;
    if (!match(r, pc + r->code[pc + 1], s, len, ofs, NULL))
      break;
  }

  *ofs = saved_offset;
}

static int is_any_of(const unsigned char *p, int len, const char *s, int *ofs) {
  int i, ch;

  ch = s[*ofs];

  for (i = 0; i < len; i++)
    if (p[i] == ch) {
      (*ofs)++;
      return 1;
    }

  return 0;
}

static int is_any_but(const unsigned char *p, int len, const char *s,
                      int *ofs) {
  int i, ch;

  ch = s[*ofs];

  for (i = 0; i < len; i++)
    if (p[i] == ch) {
      return 0;
    }

  (*ofs)++;
  return 1;
}

static int lowercase(const char *s) {
  return tolower(*(const unsigned char *)s);
}

static int casecmp(const void *p1, const void *p2, size_t len) {
  const char *s1 = (const char *)p1, *s2 = (const char *)p2;
  int diff = 0;

  if (len > 0)
    do {
      diff = lowercase(s1++) - lowercase(s2++);
    } while (diff == 0 && s1[-1] != '\0' && --len > 0);

  return diff;
}

static const char *match(const struct slre *r, int pc, const char *s, int len,
                         int *ofs, struct cap *caps) {
  int n, saved_offset;
  const char *error_string = NULL;
  int (*cmp)(const void *string1, const void *string2, size_t len0);

  while (error_string == NULL && r->code[pc] != END) {

    assert(pc < r->code_size);
    assert(pc < (int)NELEMS(r->code));

    switch (r->code[pc]) {
    case BRANCH:
      saved_offset = *ofs;
      error_string = match(r, pc + 3, s, len, ofs, caps);
      if (error_string != NULL) {
        *ofs = saved_offset;
        error_string = match(r, pc + r->code[pc + 1], s, len, ofs, caps);
      }
      pc += r->code[pc + 2];
      break;

    case EXACT:
      error_string = error_no_match;
      n = r->code[pc + 2];
      cmp = r->options & SLRE_CASE_INSENSITIVE ? casecmp : memcmp;
      if (n <= len - *ofs && !cmp(s + *ofs, r->data + r->code[pc + 1], n)) {
        (*ofs) += n;
        error_string = NULL;
      }
      pc += 3;
      break;

    case QUEST:
      error_string = NULL;
      saved_offset = *ofs;
      if (match(r, pc + 2, s, len, ofs, caps) != NULL) {
        *ofs = saved_offset;
      }
      pc += r->code[pc + 1];
      break;

    case STAR:
      error_string = NULL;
      loop_greedy(r, pc, s, len, ofs);
      pc += r->code[pc + 1];
      break;

    case STARQ:
      error_string = NULL;
      loop_non_greedy(r, pc, s, len, ofs);
      pc += r->code[pc + 1];
      break;

    case PLUS:
      if ((error_string = match(r, pc + 2, s, len, ofs, caps)) != NULL)
        break;

      loop_greedy(r, pc, s, len, ofs);
      pc += r->code[pc + 1];
      break;

    case PLUSQ:
      if ((error_string = match(r, pc + 2, s, len, ofs, caps)) != NULL)
        break;

      loop_non_greedy(r, pc, s, len, ofs);
      pc += r->code[pc + 1];
      break;

    case SPACE:
      error_string = error_no_match;
      if (*ofs < len && isspace(((unsigned char *)s)[*ofs])) {
        (*ofs)++;
        error_string = NULL;
      }
      pc++;
      break;

    case NONSPACE:
      error_string = error_no_match;
      if (*ofs < len && !isspace(((unsigned char *)s)[*ofs])) {
        (*ofs)++;
        error_string = NULL;
      }
      pc++;
      break;

    case DIGIT:
      error_string = error_no_match;
      if (*ofs < len && isdigit(((unsigned char *)s)[*ofs])) {
        (*ofs)++;
        error_string = NULL;
      }
      pc++;
      break;

    case ANY:
      error_string = error_no_match;
      if (*ofs < len) {
        (*ofs)++;
        error_string = NULL;
      }
      pc++;
      break;

    case ANYOF:
      error_string = error_no_match;
      if (*ofs < len)
        error_string =
            is_any_of(r->data + r->code[pc + 1], r->code[pc + 2], s, ofs)
                ? NULL
                : error_no_match;
      pc += 3;
      break;

    case ANYBUT:
      error_string = error_no_match;
      if (*ofs < len)
        error_string =
            is_any_but(r->data + r->code[pc + 1], r->code[pc + 2], s, ofs)
                ? NULL
                : error_no_match;
      pc += 3;
      break;

    case BOL:
      error_string = *ofs == 0 ? NULL : error_no_match;
      pc++;
      break;

    case EOL:
      error_string = *ofs == len ? NULL : error_no_match;
      pc++;
      break;

    case OPEN:
      if (caps != NULL)
        caps[r->code[pc + 1]].ptr = s + *ofs;
      pc += 2;
      break;

    case CLOSE:
      if (caps != NULL)
        caps[r->code[pc + 1]].len = (s + *ofs) - caps[r->code[pc + 1]].ptr;
      pc += 2;
      break;

    case END:
      pc++;
      break;

    default:
      fprintf(stderr, "unknown cmd (%d) at %d\n", r->code[pc], pc);
      assert(0);
      break;
    }
  }

  return error_string;
}

static const char *match2(const struct slre *r, const char *buf, int len,
                          struct cap *caps) {
  int i, ofs = 0;
  const char *error_string = error_no_match;

  if (r->anchored) {
    error_string = match(r, 0, buf, len, &ofs, caps);
  } else {
    for (i = 0; i < len && error_string != NULL; i++) {
      ofs = i;
      error_string = match(r, 0, buf, len, &ofs, caps);
    }
  }

  return error_string;
}

static const char *capture_float(const struct cap *cap, void *p, size_t len) {
  const char *fmt;
  char buf[BUF512];

  switch (len) {
  case sizeof(float):
    fmt = "f";
    break;
  case sizeof(double):
    fmt = "lf";
    break;
  default:
    return "SLRE_FLOAT: unsupported size";
  }

  sprintf(buf, "%%%d%s", cap->len, fmt);
  return sscanf(cap->ptr, buf, p) == 1 ? NULL : "SLRE_FLOAT: capture failed";
}

static const char *capture_string(const struct cap *cap, void *p, size_t len) {
  if ((int)len <= cap->len) {
    return "SLRE_STRING: buffer size too small";
  }
  memcpy(p, cap->ptr, cap->len);
  ((char *)p)[cap->len] = '\0';
  return NULL;
}

static const char *capture_int(const struct cap *cap, void *p, size_t len) {
  const char *fmt;
  char buf[BUF512];

  switch (len) {
  case sizeof(char):
    fmt = "hh";
    break;
  case sizeof(long int):
    fmt = "ld";
    break;
  default:
    return "SLRE_INT: unsupported size";
  }

  sprintf(buf, "%%%d%s", cap->len, fmt);
  return sscanf(cap->ptr, buf, p) == 1 ? NULL : "SLRE_INT: capture failed";
}

static const char *capture(const struct cap *caps, int num_caps, va_list ap) {
  int i, type;
  size_t size;
  void *p;
  const char *err = NULL;

  for (i = 0; i < num_caps; i++) {
    type = va_arg(ap, int);
    size = va_arg(ap, size_t);
    p = va_arg(ap, void *);
    switch (type) {
    case SLRE_INT:
      err = capture_int(&caps[i], p, size);
      break;
    case SLRE_FLOAT:
      err = capture_float(&caps[i], p, size);
      break;
    case SLRE_STRING:
      err = capture_string(&caps[i], p, size);
      break;
    default:
      err = "Unknown type, expected SLRE_(INT|FLOAT|STRING)";
      break;
    }
  }
  return err;
}

int lux_match(enum slre_option options, const char *re, const char *buf) {
  struct slre slre;
  struct cap caps[32];
  int buf_len = strlen(buf);

  slre.options = options;
  return compile2(&slre, re) == NULL &&
         match2(&slre, buf, buf_len, caps) == NULL;
}

int lux_ncmatch(const char *buf, const char *re) {
  return lux_match(SLRE_CASE_INSENSITIVE, re, buf);
}
