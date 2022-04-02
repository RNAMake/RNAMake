//
//  pose_factory.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 8/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/pose_factory.h"

#include <algorithm>
#include <map>

#include "motif/motif.h"
#include "structure/chain.h"
#include "util/x3dna.h"

namespace motif {

PoseOP PoseFactory::pose_from_motif_tree(
    structure::StructureOP const& structure,
    structure::BasepairOPs const& basepairs, motif::MotifOPs const& motifs,
    std::map<util::Uuid, int, util::UuidCompare> const& designable) {
  auto ends = mf_._setup_basepair_ends(structure, basepairs);
  auto m = std::make_shared<Motif>(structure, basepairs, ends);
  mf_._setup_secondary_structure(m);
  _standardize_prepose(m);

  auto p = std::make_shared<Pose>(m);
  _add_motifs_to_pose(p, motifs);
  _add_secondary_structure_motifs(p);

  return p;
}

PoseOP PoseFactory::pose_from_file(String const& path, int gu_are_helix,
                                   int signlet_bp_seperation) {
  auto base_motif = mf_.motif_from_file(path);
  auto p = std::make_shared<Pose>(base_motif);
  _setup_motifs_from_x3dna(p, gu_are_helix, signlet_bp_seperation);

  return p;
}

void PoseFactory::_setup_motifs_from_x3dna(PoseOP& p, int gu_are_helix,
                                           int singlet_bp_seperation) {
  auto pdb_path = p->path();
  auto fname = base::filename(pdb_path);
  if (base::is_dir(p->path())) {
    pdb_path = pdb_path + "/" + fname + ".pdb";
  }

  util::X3dna x;
  auto x3_motifs = x.get_motifs(pdb_path);
  motif::MotifOPs motifs;

  for (auto const& x3_m : x3_motifs) {
    auto m = _convert_x3dna_to_motif(x3_m, p);
    motifs.push_back(m);
  }

  std::map<util::MotifType, motif::MotifOPs> motif_map;
  motif_map[util::MotifType::ALL] = motif::MotifOPs();
  for (auto const& m : motifs) {
    if (motif_map.find(m->mtype()) == motif_map.end()) {
      motif_map[m->mtype()] = motif::MotifOPs();
    }
    motif_map[util::MotifType::ALL].push_back(m);
    motif_map[m->mtype()].push_back(m);
  }

  p->set_motifs(motif_map);
}

void PoseFactory::_add_secondary_structure_motifs(PoseOP& p) {
  /*auto ss = p->secondary_structure();
  int i = 0;
  String type_name;
  std::map<String, secondary_structure::motif::MotifOPs> ss_motif_map;
  ss_motif_map["ALL"] = secondary_structure::motif::MotifOPs();
  for(auto const & m : p->motifs(util::MotifType::ALL)) {
      secondary_structure::structure::BasepairOPs ss_ends, ss_bps;
      secondary_structure::structure::ChainOPs ss_chains;
      Strings ss_end_ids;
      for(auto const & c : m->chains()) {
          secondary_structure::ResidueOPs ss_res;
          for(auto const & r : c->residues()) {
              auto ss_r = ss->get_residue(r->uuid());
              ss_res.push_back(ss_r);
          }
          ss_chains.push_back(std::make_shared<secondary_structure::structure::Chain>(ss_res));
      }
      for(auto const & bp : m->basepairs()) {
          auto res1 = ss->get_residue(bp->res1()->uuid());
          auto res2 = ss->get_residue(bp->res2()->uuid());
          auto ss_bp = ss->get_basepair(res1, res2);
          if(ss_bp == nullptr) { continue; }
          ss_bps.push_back(ss_bp);
      }
      i = 0;
      for(auto const & end : m->ends()) {
          auto res1 = ss->get_residue(end->res1()->uuid());
          auto res2 = ss->get_residue(end->res2()->uuid());
          auto ss_bp = ss->get_basepair(res1, res2);
          if(ss_bp == nullptr) {
              std::runtime_error("cannot find ss end in
  PoseFactory::_add_secondary_structure_motifs");
          }
          ss_bps.push_back(ss_bp);
          ss_end_ids.push_back(m->end_ids()[i]);
          i++;
      }
      type_name = type_to_str(m->mtype());
      auto ss_motif = std::make_shared<secondary_structure::Motif>(type_name,
  ss_ends, ss_chains);
      //ss_motif->basepairs(ss_bps);
      ss_motif->name(m->name());
      ss_motif->end_ids(ss_end_ids);

      if(ss_motif_map.find(type_name) == ss_motif_map.end()) {
          ss_motif_map[type_name] = secondary_structure::motif::MotifOPs();
      }

      ss_motif_map[type_name].push_back(ss_motif);
      ss_motif_map["ALL"].push_back(ss_motif);
  }*/

  // p->set_ss_motifs(ss_motif_map);
}

void PoseFactory::_add_motifs_to_pose(PoseOP& p,
                                      motif::MotifOPs const& motifs) {
  // int j = 0;
  float dist, best_dist;
  math::Point r1_cent, r2_cent;
  structure::ResidueOP best_match;
  std::map<util::MotifType, motif::MotifOPs> motif_map;
  motif_map[util::MotifType::ALL] = motif::MotifOPs();
  for (auto const& m : motifs) {
    structure::BasepairOPs bps;
    std::map<structure::ResidueOP, structure::ResidueOP> residue_map;

    for (auto const& bp : m->basepairs()) {
      auto new_bp = p->get_basepair(bp->uuid());
      if (new_bp.size() != 0) {
        bps.push_back(new_bp[0]);
        continue;
      }

      int best = 1000;
      structure::BasepairOP best_bp = nullptr;
      for (auto const& m2 : motifs) {
        if (m == m2) {
          continue;
        }
        for (auto const& end : m2->ends()) {
          auto alt_bp = p->get_basepair(end->uuid());
          if (alt_bp.size() == 0) {
            continue;
          }
          dist = bp->d().distance(alt_bp[0]->d());
          if (dist < best) {
            best_bp = alt_bp[0];
            best = dist;
          }
        }
      }

      if (best_bp == nullptr) {
        continue;
      }
      new_bp = structure::BasepairOPs{best_bp};
      for (auto const& r1 : best_bp->residues()) {
        r1_cent = structure::center(r1->atoms());
        best_dist = 1000;
        best_match = nullptr;
        for (auto const& r2 : bp->residues()) {
          r2_cent = structure::center(r2->atoms());
          dist = r1_cent.distance(r2_cent);
          if (dist < best_dist) {
            best_dist = dist;
            best_match = r2;
          }
        }
        residue_map[best_match] = r1;
      }
      bps.push_back(new_bp[0]);
    }

    structure::ChainOPs chains;
    for (auto const& c : m->chains()) {
      structure::ResidueOPs res;
      for (auto const& r : c->residues()) {
        auto new_r = p->get_residue(r->uuid());
        if (new_r != nullptr) {
          res.push_back(new_r);
        } else if (residue_map.find(r) != residue_map.end()) {
          res.push_back(residue_map[r]);
        } else {
          throw std::runtime_error(
              "cannot find residue in PoseFactory::_add_motifs_to_pose");
        }
      }
      chains.push_back(std::make_shared<structure::Chain>(res));
    }

    if (bps.size() != m->basepairs().size()) {
      throw std::runtime_error(
          "something went horribly wrong, did not find all basepairs  in "
          "PoseFactory::_add_motifs_to_pose");
    }

    auto structure = std::make_shared<structure::Structure>(chains);
    auto ends = mf_._setup_basepair_ends(structure, bps);
    auto m_copy = std::make_shared<Motif>(structure, bps, ends);
    m_copy->mtype(m->mtype());
    m_copy->path(m->path());
    m_copy->name(m->name());
    mf_._setup_secondary_structure(m_copy);
    structure::BasepairOP best_end;
    String best_end_id;
    int i = 0;
    if (m_copy->mtype() != util::MotifType::HELIX) {
      structure::BasepairOPs best_ends;
      Strings best_end_ids;
      for (auto const& end : m->ends()) {
        best_dist = 1000;
        best_end = nullptr;
        best_end_id = "";
        for (auto const& c_end : m_copy->ends()) {
          dist = end->d().distance(c_end->d());
          if (dist < best_dist) {
            best_end = c_end;
            best_end_id = m->end_ids()[i];
            best_dist = dist;
          }
        }
        best_ends.push_back(best_end);
        best_end_ids.push_back(best_end_id);
      }
      i++;
    }

    if (motif_map.find(m_copy->mtype()) == motif_map.end()) {
      motif_map[m_copy->mtype()] = motif::MotifOPs();
    }
    motif_map[util::MotifType::ALL].push_back(m_copy);
    motif_map[m_copy->mtype()].push_back(m_copy);
  }

  p->set_motifs(motif_map);
}

void PoseFactory::_standardize_prepose(MotifOP& p) {
  p->get_beads(p->ends());
  auto added_helix = mf_.added_helix();

  for (auto& end : p->ends()) {
    auto m_added = get_aligned_motif(end, added_helix->ends()[0], added_helix);
    if (_steric_clash(p, m_added)) {
      continue;
    }
    end->flip();
  }
}

int PoseFactory::_steric_clash(MotifOP const& m1, MotifOP const& m2) {
  float dist;
  for (auto const& c1 : m1->beads()) {
    for (auto const& c2 : m2->beads()) {
      if (c1.btype() == structure::BeadType::PHOS ||
          c2.btype() == structure::BeadType::PHOS) {
        continue;
      }
      dist = c1.center().distance(c2.center());
      if (dist < clash_radius_) {
        return 1;
      }
    }
  }

  return 0;
}

MotifOP PoseFactory::_convert_x3dna_to_motif(util::X3dna::X3Motif const& xm,
                                             PoseOP const& p) {
  /*structure::ResidueOPs res;
  for (auto const & xr : xm.residues) {
      auto r = p->get_residue(xr.num, xr.chain_id, xr.i_code);
      res.push_back(r);
  }

  structure::BasepairOPs basepairs;
  for (auto const & r : res) {
      auto bps = p->get_basepair(r->uuid());
      for (auto const & bp : bps) {
          if (std::find(res.begin(), res.end(), bp->res1()) != res.end() &&
              std::find(res.begin(), res.end(), bp->res2()) != res.end() &&
              std::find(basepairs.begin(), basepairs.end(), bp) !=
  basepairs.end()) { basepairs.push_back(bp);
          }
      }
  }

  auto m = mf_.motif_from_res(res, basepairs);
  m->mtype(util::str_to_type(xm.mtype));*/
  return MotifOP(nullptr);
}

}  // namespace motif
