//
//  secondary_structure_factory.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>

//RNAMake Headers
#include "secondary_structure/secondary_structure_factory.h"
#include "secondary_structure/basepair.h"

namespace sstruct {

void
SecondaryStructureFactory::_get_basepairs(
    SS_Tree const & sstree,
    SecondaryStructureOP & ss) {
    
    BasepairOPs basepairs, ends;
    for(auto const & n : sstree) {
        if(n->data()->type() != SS_NodeData::SS_Type::SS_BP) { continue; }
        
        auto res1 = n->data()->ss_chains()[0]->residues()[0];
        auto res2 = n->data()->ss_chains()[1]->residues()[0];

        auto bp = std::make_shared<Basepair>(res1, res2, Uuid());
        basepairs.push_back(bp);
        
        if(n->parent() == nullptr) {
            ends.push_back(bp);
            continue;
        }
        if(n->children().size() == 1) {
            if(n->children()[0]->data()->type() == SS_NodeData::SS_Type::SS_SEQ_BREAK) {
                ends.push_back(bp);
                continue;
            }
        }
        if(n->parent()->data()->type() == SS_NodeData::SS_Type::SS_SEQ_BREAK) {
            ends.push_back(bp);
            continue;
        }
       
    }
    
    ss->basepairs(basepairs);
    ss->ends(ends);
}

void
SecondaryStructureFactory::_get_motifs(
    SS_Tree const & sstree,
    SecondaryStructureOP & ss) {
    
    std::map<String, MotifOPs> motifs;
    motifs["ALL"] = MotifOPs();
    
    for(auto const & n : sstree) {
        if(n->parent_index() == -1 ||
           n->parent()->data()->type() == SS_NodeData::SS_Type::SS_SEQ_BREAK) {
            continue;
        }
        
        String type_name;
        std::vector<TreeNodeOP<SS_NodeDataOP>> nodes, bp_nodes;
        if(n->data()->type() == SS_NodeData::SS_Type::SS_BP &&
           n->parent()->data()->type() == SS_NodeData::SS_Type::SS_BP) {
            nodes = {n->parent(), n};
            bp_nodes = {n->parent(), n};
            type_name = "BP_STEP";
        }
        
        else if(n->data()->type() != SS_NodeData::SS_Type::SS_SEQ_BREAK &&
                n->parent()->data()->type() == SS_NodeData::SS_Type::SS_BP) {
            nodes = {n->parent(), n};
            bp_nodes = {n->parent()};
            type_name = n->data()->what().substr(3);
            for(auto const & c : n->children()) {
                if(c->data()->type() == SS_NodeData::SS_Type::SS_BP ||
                   c->data()->type() == SS_NodeData::SS_Type::SS_PSEUDO_BP) {
                    nodes.push_back(c);
                    bp_nodes.push_back(c);
                }
                else {
                    throw std::runtime_error("unexpected connectivity in ss_tree");
                }
            }
        }
        
        else { continue; }
        
        auto chains = _get_chains(ss, nodes);
        BasepairOPs basepairs, ends;
        for(auto const & bp_n : bp_nodes) {
            auto bp = ss->get_bp(bp_n->data()->ss_chains()[0]->residues()[0],
                                 bp_n->data()->ss_chains()[1]->residues()[0]);
            ends.push_back(bp);
        }
        
        auto m = std::make_shared<Motif>(type_name, ends, chains);
        std::map<ResidueOP, int> res;
        for(auto const & r : m->residues()) { res[r] = 1; }
        for(auto const & bp : ss->basepairs()) {
            if(res.find(bp->res1()) != res.end() &&
               res.find(bp->res2()) != res.end()) {
                basepairs.push_back(bp);
            }
        }
        m->basepairs(basepairs);
        
        if(motifs.find(type_name) == motifs.end()) {
            motifs[type_name] = MotifOPs();
        }
        
        motifs[type_name].push_back(m);
        motifs["ALL"].push_back(m);
        
    }
    
    ss->set_motifs(motifs);
        
}
        
ChainOPs
SecondaryStructureFactory::_get_chains(
    SecondaryStructureOP const & ss,
    std::vector<TreeNodeOP<SS_NodeDataOP>> const & nodes) {
    
    ResidueOPs res;
    for(auto const & n : nodes) {
        for(auto const & c : n->data()->ss_chains()) {
            for(auto const & r : c->residues()) {
                res.push_back(r);
            }
        }
    }
    std::sort(res.begin(), res.end(), res_less_than_key());
    int last = -1;
    ChainOPs chains;
    ResidueOPs c_res;
    int is_chain_start = 0;
    for(auto const & r : res) {
        is_chain_start = 0;
        for(auto const & c : ss->chains()) {
            if(c->first() == r) {
                is_chain_start = 1;
                break;
            }
        }
        if(last != -1 && (last+1 != r->num() || is_chain_start)) {
            chains.push_back(std::make_shared<Chain>(c_res));
            c_res = ResidueOPs();
        }
        c_res.push_back(r);
        last = r->num();
    }
    if(c_res.size() > 0) {
        chains.push_back(std::make_shared<Chain>(c_res));
    }
    
    return chains;
    
}
    

} //sstruct
























