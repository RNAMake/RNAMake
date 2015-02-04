//
//  structure.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <fstream>
#include "structure.h"
#include "chain.h"

void
Structure::_build_chains(
    ResidueOPs & residues) {
    chains_ = Chains();
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
                chains_.push_back(Chain(current_chain_res));
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
    Chains chains;
    for (auto const & c : chains_) {
        chains.push_back(c.copy());
    }
    cstruct.chains(chains);
    cstruct._cache_coords();
    return cstruct;
}

String
Structure::to_pdb_str() {
    int acount = 1;
    String s;
    for (auto const & c : chains_) {
        s += c.to_pdb_str(acount);
    }
    return s;
}

String
Structure::to_str() {
    String s;
    for (auto const & c : chains_) {
        s += c.to_str() + ":";
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
    Chains chains;
    Strings spl = split_str_by_delimiter(s, ":");
    for( auto const & c_str : spl) {
        Chain c = str_to_chain(c_str, rts);
        chains.push_back(c);
    }
    structure.chains(chains);
    structure._cache_coords();
    return structure;
}

