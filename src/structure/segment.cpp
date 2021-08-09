
//
// Created by Joseph Yesselman on 12/20/17.
//

#include "structure/segment.h"

namespace structure {

  String
  Segment::get_pdb_str(
          int & acount,
          int & rnum,
          char & chain_id) {
      // TODO add both proteins and small molecules
      return structure_.get_pdb_str(acount, rnum, chain_id);
  }

  void
  Segment::write_pdb(
          String const & fname) const {
      std::ofstream out;
      out.open(fname.c_str());
      out << structure_.get_pdb_str(1) << std::endl;
      out.close();
  }

  void
  Segment::write_steric_beads_to_pdb(
          String const & fname) {
      structure_.write_steric_beads_to_pdb(fname);

  }

}