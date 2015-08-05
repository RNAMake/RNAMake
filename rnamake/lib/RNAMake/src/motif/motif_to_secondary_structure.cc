//
//  motif_to_secondary_structure.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include "motif/motif_to_secondary_structure.h"
#include "structure/chain.h"


sstruct::SecondaryStructureOP
MotiftoSecondaryStructure::to_secondary_structure(
    MotifOP const & motif) {
    
    BasepairOP saved_bp;
    BasepairOPs bps;
    sstruct::ChainOPs ss_chains;
    
    for(auto const & c : motif->chains()) { chains_.push_back(c); }
    open_chains_.push(chains_[0]);
    chains_.erase(chains_.begin());
    
    ChainOP best_chain;
    ResidueOP partner_r;
    String ss;
    int passes = 0;
    while (! open_chains_.empty() ) {
        auto c = open_chains_.front();
        open_chains_.pop();
        
        sstruct::ResidueOPs ss_res;
        for(auto const & r : c->residues()) {
            ss = ".";
            bps = motif->get_basepair(r->uuid());
            for(auto const & bp : bps) {
                passes = 0;
                
                if(wc_bp(bp) && bp->bp_type() == "cW-W") { passes = 1; }
                if(gu_bp(bp) && bp->bp_type() == "cW-W") { passes = 1; }
                partner_r = bp->partner(r);
                if(passes) {
                    saved_bp = bp;
                    if     (seen_bp_.find(bp) == seen_bp_.end() &&
                            seen_res_.find(r) == seen_res_.end() &&
                            seen_res_.find(partner_r) == seen_res_.end()) {
                        seen_res_[r] = 1;
                        ss = "(";
                    }
                    else if(seen_res_.find(partner_r) != seen_res_.end()) {
                        if(seen_res_[partner_r] > 1) {
                            ss = ".";
                        }
                        else {
                            ss = ")";
                            seen_res_[r] = 1;
                            seen_res_[partner_r] += 1;
                            break;
                        }
                    }
                }
                else if(seen_res_.find(r) == seen_res_.end() ) {
                    ss = ".";
                }
            }
            
            if(saved_bp != nullptr) { seen_bp_[saved_bp] = 1; }
            ss_res.push_back(std::make_shared<sstruct::Residue>(r->name(), ss, r->num(),
                                                                r->chain_id(), r->uuid(),
                                                                r->i_code()));
        }
        
        ss_chains.push_back(std::make_shared<sstruct::Chain>(ss_res));
        best_chain = _get_next_chain(motif);
        
        if(best_chain == nullptr) { break; }
        chains_.erase(std::remove(chains_.begin(), chains_.end(), best_chain), chains_.end());
        open_chains_.push(best_chain);
        
    }
    
    auto secondary_structure = std::make_shared<sstruct::SecondaryStructure>(ss_chains);
    _setup_basepairs_and_ends(secondary_structure, motif);
    
    return secondary_structure;
}


ChainOP
MotiftoSecondaryStructure::_get_next_chain(
    MotifOP const & motif) {
    
    BasepairOPs bps;
    int best_score = -1;
    int score = 0;
    for(auto const & c : chains_) {
        score = 0;
        for(auto const & r : c->residues()) {
            bps = motif->get_basepair(r->uuid());
            for(auto const & bp : bps) {
                if(seen_bp_.find(bp) != seen_bp_.end()) {
                    score += 1;
                }
            }
        }
        if(score > best_score) {
            best_score = score;
        }
    }
    
    ChainOPs best_chains;
    for(auto const & c : chains_) {
        score = 0;
        for(auto const & r : c->residues()) {
            bps = motif->get_basepair(r->uuid());
            for(auto const & bp : bps) {
                if(seen_bp_.find(bp) != seen_bp_.end()) {
                    score += 1;
                }
            }
        }
        if(score == best_score) {
            best_chains.push_back(c);
        }
    }
    
    ChainOP best_chain;
    best_score = 10000;
    int pos = 1000, i = 0;
    for(auto const & c : best_chains) {
        pos = 1000;
        i = 0;
        for(auto const & r : c->residues()) {
            bps = motif->get_basepair(r->uuid());
            for(auto const & bp : bps) {
                if(seen_bp_.find(bp) != seen_bp_.end()) {
                    pos = i;
                    break;
                }
            }
        }
        if(pos < best_score) {
            best_score = pos;
            best_chain = c;
        }
    }
    
    return best_chain;
    
}

void
MotiftoSecondaryStructure::_setup_basepairs_and_ends(
    sstruct::SecondaryStructureOP & ss,
    MotifOP const & motif) {
    
    sstruct::BasepairOPs ss_bps, ss_ends;
    for(auto const & kv : seen_bp_) {
        auto bp = kv.first;
        auto res1 = ss->get_residue(bp->res1()->uuid());
        auto res2 = ss->get_residue(bp->res2()->uuid());
        ss_bps.push_back(std::make_shared<sstruct::Basepair>(res1, res2, bp->uuid()));
    }
    ss->basepairs(ss_bps);
    
    for(auto const & end : motif->ends()) {
        auto res1 = ss->get_residue(end->res1()->uuid());
        auto res2 = ss->get_residue(end->res2()->uuid());
        auto bp = ss->get_bp(res1, res2);
        if(bp == nullptr) {
            throw std::runtime_error("did not properly find end in generating ss");
        }
        ss_ends.push_back(bp);
    }
    ss->ends(ss_ends);
    
    
}




