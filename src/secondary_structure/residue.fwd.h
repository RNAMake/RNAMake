//
// Created by Joseph Yesselman on 10/21/17.
//

#ifndef RNAMAKE_NEW_RESIDUE_FWD_H
#define RNAMAKE_NEW_RESIDUE_FWD_H

#include <vector>

namespace secondary_structure {

  class Residue;
  typedef std::shared_ptr<Residue> ResidueOP;
  typedef std::vector<ResidueOP>   ResidueOPs;
  typedef std::vector<Residue>     Residues;

}

#endif //RNAMAKE_NEW_RESIDUE_FWD_H