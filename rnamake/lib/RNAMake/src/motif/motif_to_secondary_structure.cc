//
//  motif_to_secondary_structure.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include "motif/motif_to_secondary_structure.h"
#include "structure/chain.h"


sstruct::RNAStructureOP
MotiftoSecondaryStructure::to_secondary_structure(
    RNAStructureOP const & motif) {
    
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
    
    auto struc = std::make_shared<sstruct::Structure>(ss_chains);
    return _setup_basepairs_and_ends(struc, motif);
;
}


ChainOP
MotiftoSecondaryStructure::_get_next_chain(
    RNAStructureOP const & motif) {
    
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

sstruct::RNAStructureOP
MotiftoSecondaryStructure::_setup_basepairs_and_ends(
    sstruct::StructureOP & struc,
    RNAStructureOP const & motif) {
    
    sstruct::BasepairOPs ss_bps, ss_ends;
    for(auto const & kv : seen_bp_) {
        auto bp = kv.first;
        auto res1 = struc->get_residue(bp->res1()->uuid());
        auto res2 = struc->get_residue(bp->res2()->uuid());
        ss_bps.push_back(std::make_shared<sstruct::Basepair>(res1, res2, bp->uuid()));
    }
    
    for(auto const & end : motif->ends()) {
        auto res1 = struc->get_residue(end->res1()->uuid());
        auto res2 = struc->get_residue(end->res2()->uuid());
        auto end_bp = sstruct::BasepairOP(nullptr);
        for(auto const & bp : ss_bps) {
            if(bp->res1()->uuid() == res1->uuid() &&
               bp->res2()->uuid() == res2->uuid()) {
                end_bp = bp;
                break;
            }
        }
        if(end_bp == nullptr) {
            throw std::runtime_error("did not properly find end in generating ss");
        }
        ss_ends.push_back(end_bp);
    }

    return std::make_shared<sstruct::RNAStructure>(struc, ss_bps, ss_ends);
    
}




