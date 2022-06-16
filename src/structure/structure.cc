//
// Created by Joseph Yesselman on 11/28/17.
//

#include "structure/structure.h"
#include "doctest.h"
#include "structure/pdb_parser.h"
#include <base/log.hpp>

namespace structure {

String Structure::get_pdb_str(int &acount, int &rnum, char &chain_id) const {

  auto s = String();
  auto i = -1;
  auto pos = 0;
  auto anscii_num = int(chain_id);
  for (auto const &r : _residues) {
    i++;
    s += r.get_pdb_str(acount, rnum, chain_id);
    if (_cut_points[pos] == i) {
      anscii_num += 1;
      chain_id = char(anscii_num);
    }
    rnum++;
  }
  return s;
}

void Structure::write_pdb(String const &fname) const {
  std::ofstream out;
  out.open(fname.c_str());
  out << get_pdb_str(1) << std::endl;
  out.close();
}

void Structure::write_steric_beads_to_pdb(String const &fname) {

  auto acount = 1;
  auto rnum = _residues[0].get_num();
  auto chain_id = _residues[0].get_chain_id();

  auto s = String();
  auto i = -1;
  auto pos = 0;
  auto anscii_num = int(chain_id);
  for (auto const &r : _residues) {
    i++;
    s += r.get_bead_pdb_str(acount, rnum, chain_id);
    if (_cut_points[pos] == i) {
      anscii_num += 1;
      chain_id = char(anscii_num);
    }
  }

  std::ofstream out;
  out.open(fname.c_str());
  out << s << std::endl;
  out.close();
}

// external functions
int are_residues_connected_RNA(Residue const &r1, Residue const &r2) {

  auto o3 = String("O3'");
  auto p = String("P");
  // 5' to 3'
  if (r1.get_coords(o3).distance(r2.get_coords(p)) < 5.0) {
    return 1;
  }
  // 3' to 5'
  if (r2.get_coords(o3).distance(r1.get_coords(p)) < 5.0) {
    return -1;
  }
  return 0;
}

int are_residues_connected_protein(Residue const &r1, Residue const &r2) {

  auto c = String("C");
  auto n = String("N");
  // nter to cter
  if (r1.get_coords(c).distance(r2.get_coords(n)) < 1.4) {
    return 1;
  }
  // cter to nter
  if (r2.get_coords(c).distance(r1.get_coords(n)) < 1.4) {
    return -1;
  }
  return 0;
}

int are_residues_connected(Residue const &r1, Residue const &r2) {

  if (r1.get_res_set_type() == SetType::RNA) {
    return are_residues_connected_RNA(r1, r2);
  } else if (r1.get_res_set_type() == SetType::PROTEIN) {
    return are_residues_connected_protein(r1, r2);
  } else if (r1.get_res_set_type() == SetType::UNKNOWN) {
    return 0;
  }      // small molecules are not connected into polymers
  else { // just in case I add something late
    throw StructureException("this set type is not currently supported");
  }
}

ResidueOP _get_first_residues_in_chain(ResidueOPs const &residues,
                                       std::map<ResidueOP, int> const &seen) {

  auto five_prime_end = true;
  for (auto const &r1 : residues) {
    if (seen.find(r1) != seen.end()) {
      continue;
    }
    five_prime_end = true;
    for (auto const &r2 : residues) {
      if (seen.find(r2) != seen.end()) {
        continue;
      } // not sure if I need this
      if (are_residues_connected(*r1, *r2) == -1) {
        five_prime_end = false;
        break;
      }
    }
    if (five_prime_end) {
      return r1;
    }
  }
  return nullptr;
}

ResidueOP _get_next_residue(ResidueOPs const &residues,
                            std::map<ResidueOP, int> const &seen) {
  auto lowest = ResidueOP(nullptr);
  for (auto const &r : residues) {
    if (seen.find(r) != seen.end()) {
      continue;
    }
    if (lowest == nullptr) {
      lowest = r;
    }
    if (r->get_chain_id() < lowest->get_chain_id()) {
      lowest = r;
    } else if (r->get_chain_id() == lowest->get_chain_id() &&
               r->get_num() < lowest->get_num()) {
      lowest = r;
    }
  }
  return lowest;
}

StructureOP get_structure_from_residues(ResidueOPs const &residues) {

  expects<StructureException>(
      residues.size() > 0,
      "no residues were supplied in get_structure_from_residues");

  auto current = ResidueOP(nullptr);
  auto all_res = Residues();
  auto chain_cuts = Indexes();
  auto seen = std::map<ResidueOP, int>();

  // has to be all RNA OR all Protein etc
  auto res_set_type = residues[0]->get_res_set_type();
  for (auto const &r : residues) {
    if (r->get_res_set_type() != res_set_type) {
      throw StructureException(
          "cannot generate structure with two different residue types");
    }
  }

  while (true) {
    current = _get_first_residues_in_chain(residues, seen);
    // backup get lowest residue num and chain id
    if (current == nullptr && seen.size() != residues.size()) {
      current = _get_next_residue(residues, seen);
    } else if (seen.size() == residues.size()) {
      break;
    }

    seen[current] = 1;
    auto found = true;
    while (found) {
      all_res.push_back(*current);
      found = false;
      for (auto const &r : residues) {
        if (seen.find(r) != seen.end()) {
          continue;
        }
        if (are_residues_connected(*current, *r) == 1) {
          current = r;
          found = true;
          break;
        }
      }
      if (found) {
        seen[current] = 1;
      } else {
        chain_cuts.push_back((int)all_res.size());
        break;
      }
    }
  }

  if (chain_cuts.back() != (int)all_res.size()) {
    chain_cuts.push_back((int)all_res.size());
  }

  return std::make_shared<Structure>(all_res, chain_cuts);
}

StructureOP get_structure_from_pdb(String const &pdb_name,
                                   ResidueTypeSet const &rts,
                                   SetType set_type) {

  auto pdb_parser = structure::PDBParser(rts);
  auto residues = pdb_parser.parse(pdb_name);
  if (set_type == SetType::RNA) {
    return get_structure_from_residues(residues->RNA_residues);
  } else if (set_type == SetType::PROTEIN) {
    return get_structure_from_residues(residues->protein_residues);
  } else if (set_type == SetType::UNKNOWN) {
    return get_structure_from_residues(residues->small_molecule_residues);
  } else {
    throw StructureException("residue set type not supported");
  }
}

base::VectorContainerOP<Basepair>
get_basepairs_from_x3dna(util::X3dna::X3Basepairs const &x3dna_basepairs,
                         Structure const &s) {

  auto bps = Basepairs();
  for (auto const &xbp : x3dna_basepairs) {
    Residue const *res1;
    Residue const *res2;
    try {
      res1 = &s.get_residue(xbp.res1.num, xbp.res1.chain_id, xbp.res1.i_code);
    } catch (StructureException) {

      LOGW << "cannot find RNA residue in basepair with num: " +
                  std::to_string(xbp.res1.num) +
                  " chain_id: " + String(xbp.res1.chain_id, 1) +
                  " i_code: " + String(xbp.res1.i_code, 1) +
                  "in structure SKIPPING!";
      continue;
    }

    try {
      res2 = &s.get_residue(xbp.res2.num, xbp.res2.chain_id, xbp.res2.i_code);
    } catch (StructureException) {
      LOGW << "cannot find RNA residue in basepair with num: " +
                  std::to_string(xbp.res2.num) +
                  " chain_id: " + String(xbp.res2.chain_id, 1) +
                  " i_code: " + String(xbp.res2.i_code, 1) +
                  "in structure SKIPPING!";
      continue;
    }

    // calculate center of base pair
    auto center = math::Vector3();
    auto count = 0;
    for (auto const &a : *res1) {
      center += a.get_coords();
      count += 1;
    }
    for (auto const &a : *res2) {
      center += a.get_coords();
      count += 1;
    }
    center /= count;

    auto bp_name =
        std::make_shared<base::SimpleString>(generate_bp_name(*res1, *res2));
    auto bp_type = generate_bp_type(*res1, *res2, xbp.bp_type);
    auto c1_prime_coords =
        math::Points{res1->get_coords("C1'"), res2->get_coords("C1'")};

    bps.push_back(Basepair(res1->get_uuid(), res2->get_uuid(), util::Uuid(),
                           bp_type, bp_name, xbp.bp_type, xbp.r, center,
                           c1_prime_coords));
  }
  return std::make_shared<base::VectorContainer<Basepair>>(bps);
}

} // namespace structure