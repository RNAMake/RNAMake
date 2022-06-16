//
// Created by Joseph Yesselman on 4/25/18.
//

#include "motif_data_structure/motif_state_ensemble_graph.h"

void motif_data_structure::motif_state_ensemble_graph_from_motif_graph(
    MotifGraph &mg, resources::Manager &rm,
    MotifStateEnsembleGraph &mseg /* return */,
    std::map<int, int> &index_hash /* return */) {
  int i = 0, j = 0;
  for (auto const &n : mg) {
    // build / get motif state ensemble
    auto mse = motif::MotifStateEnsembleOP(nullptr);

    if (n->data()->mtype() == util::MotifType::HELIX) {
      // is not a basepair step
      if (n->data()->residues().size() > 4) {
        LOG_ERROR << "supplied a helix motif: " + n->data()->name() +
                         " that is not a basepair "
                  << "step, this is no supported";
        exit(0);
      }

      try {
        mse = rm.motif_state_ensemble(n->data()->end_ids()[0]);
      } catch (resources::ResourceManagerException const &e) {
        LOG_ERROR << "cannot find motif state ensemble for basepair with id: " +
                         n->data()->end_ids()[0]
                  << "check to make sure its a Watson-Crick basepair step";
        exit(0);
      }
    } else {
      mse = std::make_shared<motif::MotifStateEnsemble>(n->data()->get_state());
    }

    if (i == 0) {
      j = mseg.add_ensemble(*mse);
    } else {
      int pi = index_hash[mg.parent_index(n->index())];
      int pie = mg.parent_end_index(n->index());
      j = mseg.add_ensemble(*mse, data_structure::NodeIndexandEdge{pi, pie});
    }

    index_hash[n->index()] = j;
    i++;
  }
}
