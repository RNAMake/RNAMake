//
//  structure_factory.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "structure/structure_factory.h"
#include "structure/chain.h"


StructureOP
StructureFactory::get_structure(
    String const & pdb_path) {
    
    PDBParser pdb_parser;
    auto residues = pdb_parser.parse(pdb_path);
    auto chains = build_chains(residues);
    auto s = std::make_shared<Structure>(chains);
    return s;
}

ChainOPs
StructureFactory::build_chains(
    ResidueOPs & residues) {
    ChainOPs chains;
    ResidueOP current;
    ResidueOPs current_chain_res;
    int five_prime_end = 1;
    int found = 1;
    while (true) {
        for (auto & r1 : residues) {
            five_prime_end = 1;
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
                chains.push_back(ChainOP(new Chain(current_chain_res)));
            }
        }
        if(residues.size() == 0) {
            break;
        }
    }
    
    return chains;
}