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

namespace util::x3dna {
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

void get_bp_3char_symbols(long bp_type, char zdir, char *bp_sym) {
  sprintf(bp_sym, "%c%c%c", (bp_type == 2) ? '-' : '*',
          (bp_type > 0) ? '-' : '*', zdir);
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

int case_strncmp(const char *s1, const char *s2, long n) {
  int i, c1, c2;

  for (i = 0;
       (c1 = toupper((int)s1[i]), c2 = toupper((int)s2[i]), c1 == c2) && i < n;
       i++)
    if (c1 == '\0')
      return 0;

  return (i >= n) ? 0 : c1 - c2;
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

void clear_my_globals() {
  free_cmatrix(Gvars.ATOMLIST, 1, -1, 0, -1);
  free_cmatrix(Gvars.BASELIST, 1, -1, 0, -1);
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
long number_of_atoms(const String & pdbfile, const miscPars & pars) {
  using namespace util::io;
  long n = 0;
  PDBParser parser;
  PDBLineType type;
  parser.parse(pdbfile);
  while(!parser.done()) {
    parser.next();
    type = parser.get_line_type();
    if(type == PDBLineType::ATOM) {
      n += 1;
    }
    else if(type == PDBLineType::END || type == PDBLineType::ENDMDL) {
      break;
    }
  }
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
long read_pdb(String pdbfile, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs,
              long hetatm, char *ALT_LIST) {
  char *p0, *pn;
  char rname[4], str[BUF512], str0[BUF512], temp[BUF512];
  double occupancy;
  long i, n = 0, nlen, k, occ_chk = false, modelNum = 0;
  std::ifstream in;
  String line;
  in.open(pdbfile);
  while (getline(in, line)) {
    strcpy(str, line.c_str());
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

      /*occupancy = get_occupancy(nlen, str, pdbfile);
      if (occupancy <= 0) { // ignore 0-occupancy atom 
        if (!occ_chk) {
          occ_chk = true;
          fprintf(stderr, "[i] File '%s' with atom occupancy <= 0 [%s]\n",
                  pdbfile, str);
        }
        continue;
      }    */

      n++;

      if (AtomSNum != nullptr) { /* read in original atom serial number */
        strncpy(temp, str + 6, 5);
        temp[5] = '\0'; /* 5 digits */
        if (sscanf(temp, "%5ld", &AtomSNum[n]) != 1) {
          fprintf(stderr, "Atom serial number ? %s ?\n", str);
          AtomSNum[n] = n; /* sequential number */
        }
      }

      strncpy(AtomName[n], str + 12, 4);
      AtomName[n][4] = '\0';

      if (Gvars.AtomName0 && Gvars.Name0) {
        strcpy(Gvars.AtomName0[n], AtomName[n]);
      }
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
      if (Miscs != nullptr) {
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
  if (!n)
    fprintf(stderr, "PDB file <%s> has NO ATOM/HETATM records\n", pdbfile.c_str());

  return n;
}

/* free all PDB relevant arrays to make code clean and short */
void free_pdb(long num, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs) {
  if (AtomSNum != nullptr)
    free_lvector(AtomSNum, 1, num);
  if (AtomName != nullptr)
    free_cmatrix(AtomName, 1, num, 0, 4);
  if (ResName != nullptr)
    free_cmatrix(ResName, 1, num, 0, 3);
  if (ChainID != nullptr)
    free_cvector(ChainID, 1, num);
  if (ResSeq != nullptr)
    free_lvector(ResSeq, 1, num);
  if (xyz != nullptr)
    free_dmatrix(xyz, 1, num, 1, 3);
  if (Miscs != nullptr)
    free_cmatrix(Miscs, 1, num, 0, NMISC);
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

long is_valid_base(char c, char *valid_bases) {
  if (isspace((int)c)) /* skip white space */
    return false;

  if (strchr(valid_bases, c) == NULL) {
    fprintf(stderr, "skip %c: acceptable bases [%s]\n", c, valid_bases);
    return false;
  }

  return true;
}

long repeat_num() {
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

/* assuming z1 and z2 are already normalized */
double z1_z2_angle_in_0_to_90(double *z1, double *z2) {
  double dircos = dot(z1, z2);

  return 90.0 - fabs(dot2ang(dircos) - 90.0); /* absolute value */
}

void print_frame(FILE *fp, double *O, double **R) {
  fprintf(fp, "%10.4f %10.4f %10.4f  # origin\n", O[1], O[2], O[3]);
  fprintf(fp, "%10.4f %10.4f %10.4f  # x-axis\n", R[1][1], R[2][1], R[3][1]);
  fprintf(fp, "%10.4f %10.4f %10.4f  # y-axis\n", R[1][2], R[2][2], R[3][2]);
  fprintf(fp, "%10.4f %10.4f %10.4f  # z-axis\n", R[1][3], R[2][3], R[3][3]);
}


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

int dval_compare(const void *v1, const void *v2) {
  const double *p1, *p2;

  p1 = (const double *)v1;
  p2 = (const double *)v2;

  return (*p1 > *p2) ? 1 : (*p1 < *p2) ? -1 : 0;
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

static const char *error_no_match = "No match";

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

