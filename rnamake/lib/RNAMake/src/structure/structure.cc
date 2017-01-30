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


String
Structure::to_pdb_str(
    int renumber) {
    
    int rnum = -1;
    String chain_id = "";
    
    if(renumber != -1) {
        rnum = 1;
        chain_id = "A";
    }
    
    int acount = 1;
    String s;
    for (auto const & c : chains_) {
        s += c->to_pdb_str(acount, rnum, chain_id);
        if(renumber != -1) {
            rnum += c->length();
        }
        s += "TER\n";
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
    String const fname,
    int renumber) {
    std::ofstream out;
    out.open(fname.c_str());
    String s = to_pdb_str(renumber);
    out << s << std::endl;
    out.close();
}



StructureOP
structure_from_pdb(
        String const & path,
        ResidueTypeSet const & rts) {

    auto pdb_parser = PDBParser(rts);
    auto residues = pdb_parser.parse(path);
    auto chains = ChainOPs();
    connect_residues_into_chains(residues, chains);
    return std::make_shared<Structure>(chains);
}


