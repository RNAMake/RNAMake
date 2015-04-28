//
//  mtst_to_met_converter.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 3/20/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_ensemble_tree.h"
#include "motif_ensemble.h"
#include "mtst_to_met_converter.h"
#include "motif_tree_node.h"
#include "chain.h"
#include "residue.h"


MotifEnsembleTreeOP
MTSTtoMETConverter::convert(
    MotifTreeStateTree const & mtst,
    MotifTree & mt,
    String const & seq,
    int start_pos) {
    
    mt_ = mt;
    p_ = mt.to_pose();
    
    if(seq.length() == 0) {
        dseq_ = p_->sequence();
    }
    else {
        dseq_ = seq;
    }
    met_ = MotifEnsembleTreeOP ( new MotifEnsembleTree() );
    node_num_map_.clear();
    
    ChainOP start_chain = _get_start_chain(mt_.nodes()[start_pos]);
    String motif_bp, last_bp, next_bp, last_step;
    MotifTreeNodeOP other_node, partner;
    int i = -1;
    int flip, start_index;
    MotifEnsembleTreeNodeOP parent;
    int end_index, size;
    
    //figure out the correct parent node indices since met_ will have different
    //numbering using single basepairs instead of whole helices
    int n_count = -1;
    for( auto const & n : mt_.nodes()) {
        i++;
        if(i < start_pos) {
            n_count++;
            node_num_map_[i] = n_count;
            continue;
        }
        
        if(n->motif()->mtype() == HELIX) {
            size = (int)n->motif()->residues().size();
            size /= 2;
            size -= 1;
            n_count += size;
            node_num_map_[i] = n_count;
        }
        
        else {
            n_count++;
            node_num_map_[i] = n_count;
        }
    }
    
    i = -1;
    
    for (auto const & n : mt_.nodes()) {
        i++;
        if(i == 0) { continue; }
        parent = met_->nodes()[ node_num_map_ [mtst.nodes()[i]->parent()->index() ] ] ;
        end_index = mtst.nodes()[i]->parent_end_index();
        if(i < start_pos) {
            MotifState ms (mtst.nodes()[i]->mts(), 1.0);
            MotifEnsemble me;
            me.add_motif_state(ms);
            met_->add_ensemble(me, parent, end_index);
            continue;
        }
        if(n->motif()->mtype() == HELIX) {
            if(i > start_pos) {
                motif_bp = _get_next_bp(mt_.nodes()[i-1], n, start_chain, 1);
            }
            flip = mtst.nodes()[i]->mts()->flip();
            start_index = mtst.nodes()[i]->mts()->start_index();
            last_bp = _add_helix(n, start_chain, motif_bp, start_index, flip, parent, end_index);
            
            if (i < mt_.nodes().size()-1) {
                next_bp = _get_next_bp(n, mt_.nodes()[i+1], start_chain);
                last_step = last_bp + "=" + next_bp;
                MotifEnsemble me (last_step, start_index, flip);
                met_->add_ensemble(me);
            }
            
            else if (n->connections().size() > 1) {
                other_node = NULL;
                for (auto const & c : n->connections()) {
                    partner = c->partner(n);
                    if(partner == mt_.nodes()[i-1]) { continue; }
                    other_node = partner;
                    break;
                }
                next_bp = _get_next_bp(n, other_node, start_chain);
                last_step = last_bp + "=" + next_bp;
                MotifEnsemble me (last_step, start_index, flip);
                met_->add_ensemble(me);
            }
            
        }
        
        else {
            MotifEnsemble me;
            MotifState ms (mtst.nodes()[i]->mts(), 1.0);
            me.add_motif_state(ms);
            met_->add_ensemble(me, parent, end_index);
        }
    }
    
    return met_;
}

int
MTSTtoMETConverter::_get_chain_pos(
    ResidueOP const & res) {
    int i = 0;
    for(auto const & c : p_->chains()) {
        for (auto const & r : c->residues()) {
            if(res == r) { return i; }
        }
        i++;
    }
    return -1;
}

String
MTSTtoMETConverter::_add_helix(
    MotifTreeNodeOP const & n,
    ChainOP const & chain,
    String const & motif_bp,
    int const & start_index,
    int const & flip,
    MotifEnsembleTreeNodeOP const & parent,
    int const & end_index) {
    
    Strings bps_str;
    if(motif_bp.size() > 0) { bps_str.push_back(motif_bp); }
    int chain_pos = _get_chain_pos(chain->first());
    int start_pos = 10000;
    int i = 0;
    ResidueOP r_new;
    String res1, res2;
    BasepairOPs bps;
    for (auto const & r : n->motif()->residues()) {
        r_new = p_->get_residue(r->uuid());
        if(r_new == NULL) { continue; }
        i = 0;
        for (auto const & rn : chain->residues()) {
            if ( r_new == rn && i < start_pos) { start_pos = i; }
            i++;
        }
    }
    
    i = -1;
    for(auto const & r : chain->residues()) {
        i++;
        if(i < start_pos) { continue; }
        r_new = n->motif()->get_residue(r->uuid());
        if(r_new == NULL) { break; }
        res1 = dseq_[r->num()-1 + chain_pos];
        bps = p_->get_basepair(r->uuid());
        for(auto const & bp : bps) {
            if(bp->res1() == r) {
                res2 = dseq_[bp->res2()->num()-1 + _get_chain_pos(bp->res2())];
                bps_str.push_back(res1+res2);
            }
            else {
                res2 = dseq_[bp->res1()->num()-1 + _get_chain_pos(bp->res1())];
                bps_str.push_back(res2+res1);
            }
        }
    }
    
    Strings steps;
    for(i = 1; i < bps_str.size(); i++) {
        steps.push_back(bps_str[i-1]+"="+bps_str[i]);
    }
    
    MotifEnsemble me (steps[0], start_index, flip);
    met_->add_ensemble(me, parent, end_index);
    for (i = 1; i < steps.size(); i++) {
        MotifEnsemble me (steps[i], start_index, flip);
        met_->add_ensemble(me);
    }
    
    return bps_str.back();
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

String
MTSTtoMETConverter::_get_next_bp(
    MotifTreeNodeOP const & n,
    MotifTreeNodeOP const & next_node,
    ChainOP const & chain,
    int flip) {
    
    MotifTreeConnectionOP conn;
    for (auto const & c : n->connections() ) {
        if(c->partner(n) == next_node) { conn = c; break; }
    }
    
    BasepairOP bp;
    if(!flip) { bp = conn->motif_end(next_node); }
    else      { bp = conn->motif_end(n);         }
    
    ResidueOP r_new, res1;
    int best_index = 10000;
    int i = 0;
    for (auto const & r : bp->residues()) {
        r_new = p_->get_residue(r->uuid());
        i = -1;
        for(auto const & rn : chain->residues()) {
            i++;
            if(rn == r_new && i < best_index) {
                best_index = i;
                res1 = r;
                break;
            }
        }
    }
    
    if(res1 == bp->res1()) {
        return bp->res1()->short_name() + bp->res2()->short_name();
    }
    else {
        return bp->res2()->short_name() + bp->res1()->short_name();
    }
    
}
