//
//  mtst_to_met_converter.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 3/20/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "mtst_to_met_converter.h"
#include "motif_tree_node.h"
#include "chain.h"
#include "residue.h"


MotifEnsembleTreeOP
MTSTtoMETConverter::convert(
    MotifTreeStateTree const & mtst) {
    
    mt_ = mtst.to_motiftree();
    p_ = mt_.to_pose();
    dseq_ = p_->sequence();
    ChainOP start_chain = _get_start_chain(mt_.nodes()[1]);
    
    
    return NULL;
}

ChainOP
MTSTtoMETConverter::_get_start_chain(
    MotifTreeNodeOP const & n) {
    
    ChainOP start_chain = NULL;
    ResidueOP r_new;
    ResidueOPs residues = n->motif()->residues();
    int closest_to_5prime = 1000;
    int i = 0;
    
    for(auto const & r : residues) {
        r_new = p_->get_residue(r->uuid());
        if(r_new == NULL) { continue; }
        
        for(auto const & c : p_->chains()) {
            i = 0;
            for(auto const & rn : c->residues()) {
                if(r_new == rn && i < closest_to_5prime) {
                    closest_to_5prime = i;
                    start_chain = c;
                }
                i++;
            }
        }
    }
    
    return start_chain;
    
}