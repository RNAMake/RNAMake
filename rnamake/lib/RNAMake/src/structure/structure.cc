//
//  structure.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <fstream>

//RNAMake Headers
#include "structure/structure.h"
#include "structure/residue.h"
#include "structure/chain.h"
#include "structure/chain.fwd.h"



Structure
Structure::copy() {
    ChainOPs chains;
    for (auto const & c : chains_) {
        chains.push_back(ChainOP(new Chain(c->copy())));
    }
    return Structure(chains);
}

void
Structure::renumber() {
    int i = 1, j = 0;
    auto chain_ids = split_str_by_delimiter("A B C D E F G H I J K L M", " ");
    for(auto & c : chains_) {
        for(auto & r : c->residues()) {
            r->num(i);
            r->chain_id(chain_ids[j]);
            i++;
        }
        j++;
    }
    
}

ResidueOPs const
Structure::residues() {
    ResidueOPs residues;
    for (auto & c : chains_) {
        for (auto r : c->residues()) {
            residues.push_back(r);
        }
    }
    return residues;
}

ResidueOP const
Structure::get_residue(
    int const & num,
    String const & chain_id,
    String const & i_code) {
    for( auto & c : chains_) {
        for (auto & r : c->residues() ){
            if (num == r->num() && chain_id.compare(r->chain_id()) == 0 && i_code.compare(r->i_code()) == 0) {
                return r;
            }
        }
    }
    return ResidueOP(NULL);
}

ResidueOP const
Structure::get_residue(
    Uuid const & uuid) {
    for( auto & c : chains_) {
        for (auto & r : c->residues() ){
            if ( r->uuid() == uuid) { return r; }
        }
    }
    
    return ResidueOP(NULL);
}



String
Structure::to_pdb_str() {
    int acount = 1;
    String s;
    for (auto const & c : chains_) {
        s += c->to_pdb_str(acount);
    }
    return s;
}

String
Structure::to_str() {
    String s;
    for (auto const & c : chains_) {
        s += c->to_str() + ":";
    }
    return s;
}

void
Structure::to_pdb(
    String const fname) {
    std::ofstream out;
    out.open(fname.c_str());
    String s = to_pdb_str();
    out << s << std::endl;
    out.close();
}


Structure
str_to_structure(
    String const & s,
    ResidueTypeSet const & rts) {
    ChainOPs chains;
    Strings spl = split_str_by_delimiter(s, ":");
    for( auto const & c_str : spl) {
        Chain c = str_to_chain(c_str, rts);
        chains.push_back(ChainOP(new Chain(c)));
    }
    return Structure(chains);
}

