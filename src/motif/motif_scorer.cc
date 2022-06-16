//
//  motif_scorer.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/16/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/motif_scorer.h"

namespace motif {

MotifScorer::MotifScorer()
    : bp_ref_energy_(StringFloatMap()), unpaired_pentalty_(4.0) {}

float MotifScorer::score(MotifOP const &m) {

  float score = 0;
  for (auto const &bp : m->basepairs()) {
    if (bp->bp_type().compare("cW-W") == 0) {
      score += _score_cWW_bp(bp);
    } else if (bp_ref_energy_.find(bp->bp_type()) != bp_ref_energy_.end()) {
      score += bp_ref_energy_[bp->bp_type()];
    } else {
      score += 6.1122;
    }
  }

  structure::BasepairOPs bps;
  for (auto const &r : m->residues()) {
    bps = m->get_basepair(r->uuid());
    if (bps.size() == 0) {
      score += unpaired_pentalty_;
    }
  }

  score += m->residues().size() * 0.05;

  return score;
}

float MotifScorer::_score_cWW_bp(structure::BasepairOP const &bp) {

  if (!(wc_bp(bp) || gu_bp(bp))) {
    return 2;
  }
  String bpstr = bp->res1()->short_name() + bp->res2()->short_name();
  if (bpstr.compare("GC") == 0 || bpstr.compare("CG") == 0) {
    return -2;
  } else if (bpstr.compare("AU") == 0 || bpstr.compare("UA") == 0) {
    return -1;
  } else {
    return -0.5;
  }
}

void MotifScorer::_bp_reference_energy_table() {
  bp_ref_energy_["cm-"] = 6.1122307919;
  bp_ref_energy_["cM-M"] = 6.1122307919;
  bp_ref_energy_["tW+W"] = 3.11366762268;
  bp_ref_energy_["c.+M"] = 5.68522906698;
  bp_ref_energy_[".W+W"] = 6.1122307919;
  bp_ref_energy_["tW-M"] = 2.42283130036;
  bp_ref_energy_["tm-M"] = 2.71577524698;
  bp_ref_energy_["cW+M"] = 3.33339125508;
  bp_ref_energy_[".W-W"] = 4.33166562348;
  bp_ref_energy_["cM+."] = 6.1122307919;
  bp_ref_energy_["c.-m"] = 6.1122307919;
  bp_ref_energy_["cM+W"] = 4.4042238922;
  bp_ref_energy_["tM+m"] = 6.1122307919;
  bp_ref_energy_["tM-W"] = 3.02141948251;
  bp_ref_energy_["cm-m"] = 5.12076349023;
  bp_ref_energy_["cM-W"] = 6.1122307919;
  bp_ref_energy_["cW-W"] = 0.056986982519;
  bp_ref_energy_["c.-M"] = 5.43544907015;
  bp_ref_energy_["cm+M"] = 2.7132962365;
  bp_ref_energy_["cm-M"] = 3.23361276018;
  bp_ref_energy_["...."] = 4.18066203386;
  bp_ref_energy_["cm-W"] = 4.36687710812;
  bp_ref_energy_["tM-m"] = 2.83911913314;
  bp_ref_energy_["c.-W"] = 6.1122307919;
  bp_ref_energy_["cM+m"] = 5.68522906698;
  bp_ref_energy_["cM-m"] = 3.12321871743;
}

} // namespace motif