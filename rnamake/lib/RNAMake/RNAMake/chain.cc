//
//  chain.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "chain.h"

String
Chain::to_str() const {
    String s;
    for ( auto const & r : residues_ ) {
        s += r.to_str() + ";";
    }
    return s;
}

String
Chain::to_pdb_str(int & acount) const {
    String s;
    for (auto const & r : residues_ ) {
        s += r.to_pdb_str(acount);
    }
    return s;
}

void
Chain::to_pdb(String const fname) const {
    std::ofstream out;
    out.open(fname.c_str());
    int i = 1;
    String s = to_pdb_str(i);
    out << s << std::endl;
    out.close();
}

Chain
Chain::copy() const {
    Residues residues(residues_.size());
    int i = 0;
    for (auto const & r : residues_) {
        residues[i] = r.copy();
    }
    return Chain(residues);
}


Chain
str_to_chain(
    String const & s,
    ResidueTypeSet const & rts) {
    Chain c;
    Residues residues;
    Strings spl = split_str_by_delimiter(s, ";");
    for(auto const & r_str : spl) {
        Residue r = str_to_residue(r_str, rts);
        residues.push_back(r);
    }
    c.residues(residues);
    return c;
}