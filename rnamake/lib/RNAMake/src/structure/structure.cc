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

void
Structure::_build_chains(
    ResidueOPs & residues) {
    chains_ = ChainOPs();
    ResidueOP current;
    ResidueOPs current_chain_res;
    int five_prime_end = 1;
    int found = 1;
    while (true) {
        five_prime_end = 1;
        for (auto & r1 : residues) {
            for (auto & r2 : residues) {
                if(r1->connected_to(*r2) == -1) {
                    five_prime_end = 0;
                    break;
                }
            }
            if ( five_prime_end ) { current = r1; break; }
        }
        if(!five_prime_end) { break; }
        residues.erase(std::remove(residues.begin(), residues.end(), current), residues.end());
        current_chain_res = ResidueOPs();
        found = 1;
        while ( found ) {
            current_chain_res.push_back(current);
            found = 0;
            for (auto & r : residues) {
                if ( current->connected_to(*r) == 1) {
                    current = r;
                    found = 1;
                    break;
                }
            }
            if( found ) {
                residues.erase(std::remove(residues.begin(), residues.end(), current), residues.end());
            }
            else {
                chains_.push_back(ChainOP(new Chain(current_chain_res)));
            }
        }
        if(residues.size() == 0) {
            break;
        }
    }
    
}

Structure
Structure::copy() {
    Structure cstruct;
    for (auto const & c : chains_) {
        cstruct.chains_.push_back(ChainOP(new Chain(c->copy())));
    }
    cstruct._cache_coords();
    return cstruct;
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


void
Structure::renumber() {
    int i = 1;
    for( auto & r : residues()) {
        r->num(i);
        r->chain_id("A");
        i++;
    }
    
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
    Structure structure;
    ChainOPs chains;
    Strings spl = split_str_by_delimiter(s, ":");
    for( auto const & c_str : spl) {
        Chain c = str_to_chain(c_str, rts);
        chains.push_back(ChainOP(new Chain(c)));
    }
    structure.chains(chains);
    structure._cache_coords();
    return structure;
}

