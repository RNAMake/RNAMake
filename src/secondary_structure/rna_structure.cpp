//
//  rna_structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/rna_structure.h"
#include <cassert>

namespace secondary_structure {

BasepairOPs RNAStructure::get_basepair(util::Uuid const &bp_uuid) {
  BasepairOPs bps;
  for (auto const &bp : basepairs_) {
    if (bp->uuid() == bp_uuid) {
      bps.push_back(bp);
    }
    if (bp->res1()->uuid() == bp_uuid || bp->res2()->uuid() == bp_uuid) {
      bps.push_back(bp);
    }
  }
  return bps;
}

BasepairOPs RNAStructure::get_basepair(ResidueOP const &res1,
                                       ResidueOP const &res2) {
  BasepairOPs bps;
  for (auto &bp : basepairs_) {
    if (bp->res1()->uuid() == res1->uuid() &&
        bp->res2()->uuid() == res2->uuid()) {
      bps.push_back(bp);
    }
    if (bp->res1()->uuid() == res2->uuid() &&
        bp->res2()->uuid() == res1->uuid()) {
      bps.push_back(bp);
    }
  }
  return bps;
}

BasepairOPs RNAStructure::get_basepair(util::Uuid const &uuid1,
                                       util::Uuid const &uuid2) {

  BasepairOPs bps;
  for (auto const &bp : basepairs_) {
    if (bp->res1()->uuid() == uuid1 && bp->res2()->uuid() == uuid2) {
      bps.push_back(bp);
    }
    if (bp->res1()->uuid() == uuid2 && bp->res2()->uuid() == uuid1) {
      bps.push_back(bp);
    }
  }
  return bps;
}

BasepairOPs RNAStructure::get_basepair(String const &name) {
  for (auto const &bp : basepairs_) {
    if (name.compare(bp->name()) == 0) {
      return BasepairOPs{bp};
    }
  }

  throw std::runtime_error("could not find basepair with name " + name);
}

BasepairOP RNAStructure::get_end(String const &name) {
  for (auto const &bp : ends_) {
    if (name.compare(bp->name()) == 0) {
      return bp;
    }
  }

  throw std::runtime_error("could not find basepair with name " + name);
}

void RNAStructure::replace_sequence(String const &seq) {
  auto spl = base::split_str_by_delimiter(seq, "&");
  auto seq2 = String();
  for (auto const &s : spl) {
    seq2 += s;
  }

  if (spl.size() != chains().size()) {
    throw Exception("cannot replace sequence with one with a differnt number "
                    "of chains: \n org: " +
                    sequence() + "\n new: " + seq);
  }

  if (seq2.length() != residues().size()) {
    throw Exception(
        "cannot replace sequence with a different length sequence: \n org: " +
        sequence() + "\n new: " + seq);
  }

  // assert(seq2.length() == residues().size() && "cannot replace sequence with
  // a different length sequence");
  int i = 0;
  for (auto &r : residues()) {
    r->name(String(1, seq2[i]));
    i++;
  }
}

} // namespace secondary_structure
