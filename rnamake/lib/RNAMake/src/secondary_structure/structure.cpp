//
//  structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/structure.h"

namespace sstruct {

void
Structure::_setup_chains(
    String const & sequence,
    String const & dot_bracket) {
    
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
    int count = 1;
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
            chains_.push_back(std::make_shared<Chain>(res));
            res = ResidueOPs();
        }
    }
    
    if(res.size() > 0) {
        chains_.push_back(std::make_shared<Chain>(res));
    }
    
}
    

ResidueOP
Structure::get_residue(
    int const & num,
    String const & chain_id,
    String const & i_code) {
    for( auto & c : chains_) {
        for (auto & r : c->residues() ){
            if (num == r->num() &&
                chain_id.compare(r->chain_id()) == 0 &&
                i_code.compare(r->i_code()) == 0) {
                return r;
            }
        }
    }
    return ResidueOP(NULL);
}

ResidueOP 
Structure::get_residue(
    Uuid const & uuid) {
    for( auto & c : chains_) {
        for (auto & r : c->residues() ){
            if ( r->uuid() == uuid) { return r; }
        }
    }
    
    return ResidueOP(NULL);
}
 
}
