//
//  motif_factory.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/motif_factory.h"
#include "structure/chain.h"
#include "util/file_io.h"
#include "util/x3dna.h"


MotifOP
MotifFactory::motif_from_file(
    String const & path) {
    
    auto fname = filename(path);
    StructureOP structure;
    if(is_dir(path)) {
        structure = sf_.get_structure(path + "/" + fname + ".pdb");
    }
    else {
        structure = sf_.get_structure(path);
        fname = fname.substr(0, -4);
    }
    
    auto basepairs = _setup_basepairs(path, structure);
    auto ends = _setup_basepair_ends(structure, basepairs);
    auto m = std::make_shared<Motif>(structure, basepairs, ends);
    
    
    return m;
}


BasepairOPs
MotifFactory::_setup_basepairs(
    String const & path,
    StructureOP const & structure) {
    
    BasepairOPs basepairs;
    X3dna x3dna_parser;
    X3Basepairs x_basepairs = x3dna_parser.get_basepairs(path);
    ResidueOP res1, res2;
    BasepairOP bp;
    for(auto const & xbp : x_basepairs) {
        res1 = structure->get_residue(xbp.res1.num, xbp.res1.chain_id, xbp.res1.i_code);
        res2 = structure->get_residue(xbp.res2.num, xbp.res2.chain_id, xbp.res2.i_code);
        if (res1 == nullptr || res2 == nullptr) {
            throw "cannot find residues in basepair during setup";
        }
        
        bp = BasepairOP(new Basepair(res1, res2, xbp.r, xbp.bp_type));
        basepairs.push_back(bp);
    }
    
    return basepairs;
    
}


BasepairOPs
MotifFactory::_setup_basepair_ends(
    StructureOP const & structure,
    BasepairOPs const & basepairs) {
    
    ResidueOPs chain_ends;
    for(auto const & c : structure->chains()) {
        chain_ends.push_back(c->first());
        if(c->residues().size() > 1) { chain_ends.push_back(c->last()); }
    }
    
    BasepairOPs ends;
    for( auto const & bp : basepairs ) {
        for (auto const & ce1 : chain_ends) {
            for(auto const & ce2 : chain_ends) {
                if(bp->bp_type().compare("cW-W") == 0 && bp->res1() == ce1 && bp->res2() == ce2) {
                    ends.push_back(bp);
                }
            }
        }
    }
    
    return ends;
}

























