//
//  util.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <queue>

//RNAMake Headers
#include "secondary_structure/util.h"

namespace sstruct {


String
assign_end_id(
    RNAStructureOP const & ss,
    BasepairOP const & end) {
    
    
    int found = 0;
    for(auto const & e : ss->ends()) {
        if(end->uuid() == e->uuid()) { found = 1; break; }
    }
    
    if(!found) {
        throw SecondaryStructureException("supplied an end that is not in current ss element");
    }
    
    ChainOPs all_chains;
    for(auto const & c : ss->chains()) { all_chains.push_back(c); }
    std::queue<ChainOP> open_chains;
    
    for(auto const & c : all_chains) {
        if(c->first() == end->res1() || c->first() == end->res2()) {
            open_chains.push(c);
            break;
        }
    }
    
    all_chains.erase(std::remove(all_chains.begin(), all_chains.end(), open_chains.front()),
                     all_chains.end());
    
    std::map<ResidueOP, int> seen_res;
    std::map<BasepairOP, int> seen_bp;
    BasepairOP saved_bp, bp;
    BasepairOPs bps;
    String structure, sequence, dot_bracket;
    using IdChain = std::pair<String, String>;
    std::vector<IdChain> id_chains;
    int score = 0, pos = 0, i =0;
    ChainOPs best_chains;
    ChainOP best_chain;
    while( ! open_chains.empty() ) {
        ChainOP c = open_chains.front();
        open_chains.pop();
        
        for(auto const & r : c->residues()) {
            dot_bracket = ".";
            bp = nullptr;
            bps = ss->get_basepair(r->uuid());
            if(bps.size() > 0) { bp = bps[0]; }
        
            saved_bp = nullptr;
            if(bp != nullptr) {
                saved_bp = bp;
                auto partner_r = bp->partner(r);
                if     (seen_bp.find(bp) == seen_bp.end() &&
                        seen_res.find(r) == seen_res.end() &&
                        seen_res.find(partner_r) == seen_res.end()) {
                    seen_res[r] = 1;
                    dot_bracket = "(";
                }
                else if(seen_res.find(partner_r) != seen_res.end()) {
                    if(seen_res[r] > 1) {
                        dot_bracket = ".";
                    }
                    else {
                        dot_bracket = ")";
                        seen_res[r] = 1;
                        seen_res[partner_r] += 1;
                    }
                }
            }
            structure += dot_bracket;
            sequence  += r->name();
            
            if(saved_bp != nullptr) { seen_bp[saved_bp] = 1; }
        }
        
        id_chains.push_back(IdChain(sequence, structure));
        sequence = ""; structure = "";
        int best_score = -1;
        
        for(auto const & c : all_chains) {
            score = 0;
            for(auto const & r : c->residues()) {
                bps = ss->get_basepair(r->uuid());
                if(bps.size() > 0 && seen_bp.find(bps[0]) != seen_bp.end()) {
                    score += 1;
                }
            }
            if(score > best_score) { best_score = score; }
        }
        
        best_chains.resize(0);
        for(auto const & c : all_chains) {
            score = 0;
            for(auto const & r : c->residues()) {
                bps = ss->get_basepair(r->uuid());
                if(bps.size() > 0 && seen_bp.find(bps[0]) != seen_bp.end()) {
                    score += 1;
                }
            }
            if (score == best_score) {
                best_chains.push_back(c);
            }
        }
        
        best_score = 10000;
        best_chain = nullptr;
        for(auto const & c : best_chains) {
            pos = 1000;
            i = 0;
            for(auto const & r : c->residues()) {
                bps = ss->get_basepair(r->uuid());
                if(bps.size() > 0 && seen_bp.find(bps[0]) != seen_bp.end()) {
                    pos = i;
                    break;
                }
                i++;
            }
            if(pos < best_score) {
                best_score = pos;
                best_chain = c;
            }
        }
        
        if(best_chain == nullptr) { break; }
        all_chains.erase(std::remove(all_chains.begin(), all_chains.end(), best_chain),
                         all_chains.end());
        open_chains.push(best_chain);
    }
    
    String ss_id = "";
    i = 0;
    for(auto const & id_chain : id_chains) {
        ss_id += id_chain.first + "_";
        for(auto const & e : id_chain.second) {
            if     (e == '(') {
                ss_id += "L";
            }
            else if(e == ')') {
                ss_id += "R";
            }
            else if(e == '.') {
                ss_id += "U";
            }
            else {
                throw SecondaryStructureException("unexpected symbol in dot bracket notation: " + std::to_string(e));
            }
        }
        if(i != id_chains.size()-1) { ss_id += "_"; }
        i++;
    }
    
    return ss_id;
    
}
    
}
