//
//  motif_to_secondary_structure.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include "motif/motif_to_secondary_structure.h"
#include "structure/chain.h"


secondary_structure::RNAStructureOP
MotiftoSecondaryStructure::to_secondary_structure(
    structure::RNAStructureOP const & motif) {
    
    structure::BasepairOP saved_bp;
    structure::BasepairOPs bps;
    secondary_structure::ChainOPs ss_chains;
    
    reset();
    for(auto const & c : motif->chains()) { chains_.push_back(c); }
    open_chains_.push(chains_[0]);
    chains_.erase(chains_.begin());
    
    structure::ChainOP best_chain;
    structure::ResidueOP partner_r;
    String ss;
    int passes = 0;
    while (! open_chains_.empty() ) {
        auto c = open_chains_.front();
        open_chains_.pop();
        
        secondary_structure::ResidueOPs ss_res;
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
                    if     (seen_bp_.find(bp->uuid()) == seen_bp_.end() &&
                            seen_res_.find(r->uuid()) == seen_res_.end() &&
                            seen_res_.find(partner_r->uuid()) == seen_res_.end()) {
                        seen_res_[r->uuid()] = 1;
                        ss = "(";
                    }
                    else if(seen_res_.find(partner_r->uuid()) != seen_res_.end()) {
                        if(seen_res_[partner_r->uuid()] > 1) {
                            ss = ".";
                        }
                        else {
                            ss = ")";
                            seen_res_[r->uuid()] = 1;
                            seen_res_[partner_r->uuid()] += 1;
                            break;
                        }
                    }
                }
                else if(seen_res_.find(r->uuid()) == seen_res_.end() ) {
                    ss = ".";
                }
            }
            
            if(saved_bp != nullptr) { seen_bp_[saved_bp->uuid()] = saved_bp; }
            ss_res.push_back(std::make_shared<secondary_structure::Residue>(r->short_name(), ss, r->num(),
                                                                r->chain_id(), r->uuid(),
                                                                r->i_code()));
        }
        
        ss_chains.push_back(std::make_shared<secondary_structure::Chain>(ss_res));
        best_chain = _get_next_chain(motif);
        
        if(best_chain == nullptr) { break; }
        chains_.erase(std::remove(chains_.begin(), chains_.end(), best_chain), chains_.end());
        open_chains_.push(best_chain);
        
    }
    
    auto struc = std::make_shared<secondary_structure::Structure>(ss_chains);
    for(auto const & r : struc->residues()) {
    }
    
    return _setup_basepairs_and_ends(struc, motif);

}

structure::ChainOP
MotiftoSecondaryStructure::_get_next_chain(
    structure::RNAStructureOP const & motif) {
    
    structure::BasepairOPs bps;
    int best_score = -1;
    int score = 0;
    for(auto const & c : chains_) {
        score = 0;
        for(auto const & r : c->residues()) {
            bps = motif->get_basepair(r->uuid());
            for(auto const & bp : bps) {
                if(seen_bp_.find(bp->uuid()) != seen_bp_.end()) {
                    score += 1;
                }
            }
        }
        if(score > best_score) {
            best_score = score;
        }
    }
    
    structure::ChainOPs best_chains;
    for(auto const & c : chains_) {
        score = 0;
        for(auto const & r : c->residues()) {
            bps = motif->get_basepair(r->uuid());
            for(auto const & bp : bps) {
                if(seen_bp_.find(bp->uuid()) != seen_bp_.end()) {
                    score += 1;
                }
            }
        }
        if(score == best_score) {
            best_chains.push_back(c);
        }
    }
    
    structure::ChainOP best_chain;
    best_score = 10000;
    for(auto const & c : best_chains) {
        int pos = 1000;
        int i = 0;
        for(auto const & r : c->residues()) {
            bps = motif->get_basepair(r->uuid());
            for(auto const & bp : bps) {
                if(seen_bp_.find(bp->uuid()) != seen_bp_.end()) {
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

secondary_structure::RNAStructureOP
MotiftoSecondaryStructure::_setup_basepairs_and_ends(
    secondary_structure::StructureOP & struc,
    structure::RNAStructureOP const & motif) {
    
    secondary_structure::BasepairOPs ss_bps, ss_ends;
    for(auto const & kv : seen_bp_) {
        auto bp = kv.second;
        auto res1 = struc->get_residue(bp->res1()->uuid());
        auto res2 = struc->get_residue(bp->res2()->uuid());
        
        if(res1 == nullptr || res2 == nullptr) {
            throw std::runtime_error("did not properly find residues for basepairs in ss");
        }
        ss_bps.push_back(std::make_shared<secondary_structure::Basepair>(res1, res2, bp->uuid()));
    }
    
    for(auto const & end : motif->ends()) {
        auto res1 = struc->get_residue(end->res1()->uuid());
        auto res2 = struc->get_residue(end->res2()->uuid());
        auto end_bp = secondary_structure::BasepairOP(nullptr);
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

    return std::make_shared<secondary_structure::RNAStructure>(struc, ss_bps, ss_ends);
    
}




