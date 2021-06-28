//
//  motif_ensemble.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/file_io.h"
#include "base/log.h"
#include "util/csv.h"
#include "motif/motif_ensemble.h"
#include "motif/motif_factory.h"

namespace motif {

MotifStateEnsembleOP
MotifEnsemble::get_state() {
  auto motif_states = MotifStateOPs();
  auto energies = Floats();

  for(auto const & mem : members_) {
    motif_states.push_back(mem->motif->get_state());
    energies.push_back(mem->energy);
  }

  return std::make_shared<MotifStateEnsemble>(motif_states, energies);

}

void
motif_ensemble_from_csv_file(
  String const & csv_file,
  std::vector<MotifEnsembleOP> & mes) {
  auto mf = motif::MotifFactory();
  auto in = io::CSVReader<4>(csv_file);
  String path, end_0, end_1, end_2;
  in.read_header(io::ignore_missing_column, "path", "end_0", "end_1", "end_2");
  auto end_0_motifs = MotifOPs();
  auto end_1_motifs = MotifOPs();
  auto end_2_motifs = MotifOPs();
  auto weights = Floats();
  LOG_INFO << "in ensemble file: " << csv_file;
  while(in.read_row(path, end_0, end_1, end_2)) {
    if(!base::file_exists(path)) {
      LOG_ERROR << "invalid pdb path: " << path;
      exit(1);
    }
    LOG_INFO << path << " " << end_0 << " "  << end_1;
    auto motifs = motif::get_standardize_motifs(mf, path);
    auto found_0 = false, found_1 = false, found_2 = false;
    for(auto const & m : motifs) {
      if(m->end_name(0) == end_0) {
        end_0_motifs.push_back(m);
        found_0 = true;
      }
      if(m->end_name(0) == end_1) {
        end_1_motifs.push_back(m);
        found_1 = true;
      }
      if(m->end_name(0) == end_2) {
        end_2_motifs.push_back(m);
        found_2 = true;
      }
    }
    if(!found_0 && !end_0.empty()) {
      LOG_ERROR << "cannot find end " << end_0 << " in motif " << path;
      exit(1);
    }
    if(!found_1 && !end_1.empty()) {
      LOG_ERROR << "cannot find end " << end_1 << " in motif " << path;
      exit(1);
    }
    if(!found_2 && !end_2.empty()) {
      LOG_ERROR << "cannot find end " << end_2 << " in motif " << path;
      exit(1);
    }
    weights.push_back(1);
  }
  mes.push_back(std::make_shared<motif::MotifEnsemble>(end_0_motifs, weights));
  mes.push_back(std::make_shared<motif::MotifEnsemble>(end_1_motifs, weights));
  if(!end_2_motifs.empty()) {
    mes.push_back(std::make_shared<motif::MotifEnsemble>(end_2_motifs, weights));
  }

}

}