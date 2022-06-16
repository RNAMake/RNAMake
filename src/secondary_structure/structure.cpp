
//
//  structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <secondary_structure/structure.h>

namespace secondary_structure {

StructureOP get_structure_from_secondary_structure(String const &sequence,
                                                   String const &dot_bracket) {

  expects<StructureException>(
      sequence.length() == dot_bracket.length() && sequence.length() > 0,
      "cannot construct new SecondaryStructure object: new sequence and dot "
      "bracket are not the same length");

  auto residues = Residues();
  auto cut_points = Cutpoints();
  auto res_num = 1;
  auto ci = 0;
  auto chain_ids = String("ABCDEFGHIJKLMNOPQRSTUVWXZ");
  auto valid_seq = String("AGUCTN&+-");
  auto valid_ss = String("[{(.)}]");
  auto res_insertion_code = ' ';

  for (auto i = 0; i < sequence.size(); i++) {
    if (sequence[i] != '&' && sequence[i] != '+') {
      expects<StructureException>(
          base::is_char_in_string(sequence[i], valid_seq) &&
              base::is_char_in_string(dot_bracket[i], valid_ss),
          "encountered invalid sequence or secondary structure element");

      residues.push_back(Residue(sequence[i], dot_bracket[i], res_num,
                                 chain_ids[ci], res_insertion_code,
                                 util::Uuid()));
      res_num += 1;
    } else {
      expects<StructureException>(sequence[i] == dot_bracket[i],
                                  "chain breaks must occur at the same spot in "
                                  "sequence and dot_bracket");

      cut_points.push_back((int)residues.size());
      ci += 1;
      if (ci == chain_ids.length() - 1) {
        ci = 0;
      }
      res_num = 1;
    }
  }
  cut_points.push_back((int)residues.size());
  return std::make_shared<Structure>(residues, cut_points);
}

} // namespace secondary_structure