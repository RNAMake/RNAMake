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
    for (auto const & c : get_chains()) {
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
    for (auto const & r : residues_) {
        s += r->to_str() + ";";
    }
    for (auto const & i : chain_cuts_) {
        s += std::to_string(i) + " ";
    }
    s += ";";
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

bool
are_structures_equal(
        StructureOP const & s1,
        StructureOP const & s2,
        int check_uuids) {

    if (s1->num_chains() != s2->num_chains()) { return false; }

    for (int i = 0; i < s1->num_residues(); i++) {
        auto result = are_residues_equal(s1->get_residue(i), s2->get_residue(i), check_uuids);
        if (!result) { return false; }
    }

    return true;

}

StructureOP
structure_from_pdb(
        String const & path,
        ResidueTypeSet const & rts) {

    auto pdb_parser = PDBParser(rts);
    auto residues = pdb_parser.parse(path);
    auto residues_and_chain_cuts = get_chain_cuts(residues);
    return std::make_shared<Structure>(residues_and_chain_cuts.residues,
                                       residues_and_chain_cuts.chain_cuts);
}


