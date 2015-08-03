//
//  secondary_structure.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <memory.h>

#include "secondary_structure/chain.h"
#include "secondary_structure/secondary_structure.h"

namespace sstruct {


SecondaryStructure::SecondaryStructure(
    String const & sequence,
    String const & dot_bracket):
    Motif()
    {
    
    if(sequence.length() != dot_bracket.length()) {
        throw SecondaryStructureException("cannot construct new SecondaryStructure object: new sequence and dot bracket are not the same length");
    }

    if(sequence.length() == 0) {
        throw SecondaryStructureException("cannot construct new SecondaryStructure object: sequence is of lenght zero!");
    }
    
    if(dot_bracket[0] != '(' && dot_bracket[0] != '.' &&
       dot_bracket[0] != '&' && dot_bracket[0] != '+') {
        throw SecondaryStructureException("cannot construct new SecondaryStructure object: dot bracket notation for secondary structure is not valid. perhaps you flipped seq and ss?");
    }
    
    ResidueOPs res;
    int count = 0;
    String chain_ids = "ABCDEFGHIKLMNO";
    int i = -1, ci = 0;
    for(auto & s : sequence) {
        i++;
        if(s != '&' && s != '+') {
            String name = "", db = "", chain_id = "";
            name += s; db += dot_bracket[i]; chain_id += chain_ids[ci];
            auto r = std::make_shared<Residue>(name, db, count, chain_id, Uuid());
            res.push_back(r);
            count++;
        }
        else {
            ci += 1;
            this->chains_.push_back(std::make_shared<Chain>(res));
            res = ResidueOPs();
        }
    }
    
    if(res.size() > 0) {
        this->chains_.push_back(std::make_shared<Chain>(res));
    }
    
}

SecondaryStructure::SecondaryStructure(
    ChainOPs const & chains) {
 
    this->chains_ = chains;

}

    
String
assign_end_id(
    SecondaryStructureOP const &,
    BasepairOP const &) {
    
    
    
    return "";
    
}

    
}