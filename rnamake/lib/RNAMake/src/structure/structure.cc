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
Structure::renumber() {
    int i = 1, j = 0;
    auto chain_ids = base::split_str_by_delimiter("A B C D E F G H I J K L M", " ");
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
Structure::residues() const {
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
Structure::to_pdb_str(
    int renumber,
    int conect_statements) {
    
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
            rnum += (int)c->residues().size();
        }
        s += "TER\n";
    }

    if(conect_statements) {
        int a_count = 0;
        int o3_pos = 9;
        int p_pos = 1;
        for(auto const & c : chains_) {
            for(int i = 0; i < c->residues().size()-1; i++) {
                int r_size = c->residues()[i]->atoms().size();
                if(!c->residues()[i]->connected_to(*c->residues()[i+1], 1.8)) {
                    s += "CONECT " + std::to_string(a_count + o3_pos) + " " + std::to_string(a_count + r_size + p_pos) +
                         "\n";
                }
                a_count += r_size;
            }
            a_count += c->last()->atoms().size();
        }

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
        int renumber,
        int conect_statements) {
    std::ofstream out;
    out.open(fname.c_str());
    String s = to_pdb_str(renumber, conect_statements);
    out << s << std::endl;
    out.close();
}


