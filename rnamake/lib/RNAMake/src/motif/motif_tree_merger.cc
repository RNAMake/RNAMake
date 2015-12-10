//
//  motif_tree_merger.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "structure/chain.h"
#include "structure/structure.h"
#include "motif/motif_tree_merger.h"

PoseOP
MotifTreeMerger::merge(GraphStatic<MotifOP> const & graph) {
    
    if(graph.size() == 1) {
        auto m = std::make_shared<Motif>(*graph.get_node(0)->data());
        return std::make_shared<Pose>(m);
    }

    for(auto const & n : graph)  { nodes_.push_back(n); }
    
    int i = -1;
    for( auto const & n : nodes_) {
        i++;
        for( auto const & c : n->data()->chains()) {
            chains_.push_back(c->subchain(0, (int)c->length()));
        }
    }
    
    MotifTreeNodeOP start_node = nodes_[0];
    _merge_chains_in_node(start_node);
    PoseOP new_pose = _build_pose();
    
    seen_connections_.clear();
    chains_.resize(0);
    nodes_.resize(0);
    
    return new_pose;
}

PoseOP
MotifTreeMerger::_build_pose() {
    ChainOPs new_chains;
        
    for( auto const & c : chains_) {
        new_chains.push_back( ChainOP( new Chain(*c)) );
    }
    
    auto new_structure = std::make_shared<Structure>(new_chains);
    new_structure->renumber();
    
    auto residues = new_structure->residues();
    std::map<Uuid, ResidueOP, UuidCompare> uuids;
    for(auto const & res : residues) { uuids[res->uuid() ] = res; }
    
    BasepairOPs basepairs;
    MotifOPs motifs;
    std::map<Uuid, int, UuidCompare> designable;
    for(auto const & node : nodes_) {
        motifs.push_back(node->data());
        for(auto const & bp : node->data()->basepairs()) {
            if(uuids.find(bp->res1()->uuid()) ==  uuids.end()) { continue; }
            if(uuids.find(bp->res2()->uuid()) ==  uuids.end()) { continue; }
            BasepairOP cbp (new Basepair(bp->copy()));
            cbp->res1( uuids[bp->res1()->uuid()]);
            cbp->res2( uuids[bp->res2()->uuid()]);
            if(node->data()->mtype() == HELIX) { designable[cbp->uuid()] = 1; }
            basepairs.push_back(cbp);
        }
    }

    PoseOP p = pf_.pose_from_motif_tree(new_structure, basepairs, motifs, designable);
    return p;
    
}

void
MotifTreeMerger::_merge_chains_in_node(
    MotifTreeNodeOP const & node) {
    for(auto const & c : node->connections() ) {
        if(c == nullptr) { continue; }
        
        if(seen_connections_.find(c) != seen_connections_.end() ) { continue; }
        seen_connections_[c] = 1;
        auto partner = c->partner(node->index());
        auto node_chains    = _get_chains_from_connection(node, c);
        auto partner_chains = _get_chains_from_connection(partner, c);


        ChainOPs merged_chains;
        if     (node->data()->mtype() != MotifType::HELIX &&
                partner->data()->mtype() == MotifType::HELIX) {
            merged_chains = _helix_merge(node_chains, partner_chains);
        }
        else if(node->data()->mtype() == MotifType::HELIX &&
                partner->data()->mtype() != MotifType::HELIX) {
            merged_chains = _helix_merge(partner_chains, node_chains);
        }
        else {
            merged_chains = _non_helix_merge(node_chains, partner_chains);
        }
        
        std::map<ChainOP, int> used_chains;
        ChainOPs new_chains;
        for (auto const & chain : node_chains.chains()) { used_chains[chain] = 1; }
        for (auto const & chain : partner_chains.chains()) { used_chains[chain] = 1; }
        for (auto const & chain : chains_) {
            if (used_chains.find(chain) == used_chains.end()) { new_chains.push_back(chain); }
        }
        for (auto const & chain : merged_chains) {
            if (chain.get() != NULL) { new_chains.push_back(chain); }
        }
        chains_ = new_chains;
        _merge_chains_in_node(partner);
    }
    
}

ChainEndPairMap
MotifTreeMerger::_get_chains_from_connection(
    MotifTreeNodeOP const & node,
    MotifTreeConnectionOP const & c) {
    
    auto end_index = c->end_index(node->index());
    auto end = node->data()->ends()[end_index];
    return _find_chains_for_end(end);
    
}

ChainEndPairMap
MotifTreeMerger::_find_chains_for_end(
    BasepairOP const & end) {
    
    std::vector<ChainInfo> chain_infos;
    std::map<ResidueOP, Ints> seen;
    int ci_index = 0;
    int i = -1;
    
    
    for( auto const & c : chains_) {
        i++;
        for( auto const & res : end->residues() ) {
            if(seen.find(res) == seen.end() ) { seen[res] = Ints(); }
            if     (res == c->first() ) {
                chain_infos.push_back(ChainInfo(c, 0, i));
                seen[res].push_back(ci_index);
                ci_index += 1;
            }
            else if(res == c->last() ) {
                chain_infos.push_back(ChainInfo(c, 1, i));
                seen[res].push_back(ci_index);
                ci_index += 1;
            }
        }
    }
    
    if(chain_infos.size() != 2) {
        std::cout << chain_infos.size() << " " << chains_.size() << std::endl;
        throw "Could not find chain for end";
    }
    
    if(chain_infos[0].pos == 0) {
        return ChainEndPairMap(chain_infos[0].chain, chain_infos[1].chain);
    }
    else {
        return ChainEndPairMap(chain_infos[1].chain, chain_infos[0].chain);
    }
    
}

ChainOPs
MotifTreeMerger::_helix_merge(
    ChainEndPairMap const & nc,
    ChainEndPairMap const & pc) {
    
    ChainOPs merged_chains(2);
    if     ( nc.is_hairpin() && pc.is_hairpin() ) {
        merged_chains[0] = _get_merged_chain(nc.p5_chain, pc.p3_chain, 1, 1);
        //throw "cannot merge an hairpin with another hairpin"
    }
    else if( nc.is_hairpin() ) {
        ChainOP p3_chain = pc.p3_chain->subchain(0, pc.p5_chain->length()-1);
        ChainOP p5_chain = pc.p5_chain->subchain(1, pc.p5_chain->length());
        merged_chains[0] = _get_merged_hairpin(p3_chain, p5_chain, nc.p5_chain);
    }
    else if( pc.is_hairpin() ) {
        ChainOP p5_chain = pc.p5_chain->subchain(1, pc.p5_chain->length()-1);
        merged_chains[0] = _get_merged_hairpin(nc.p5_chain, nc.p3_chain, p5_chain, 1);
        
    }
    else {
        merged_chains[0] = _get_merged_chain(nc.p5_chain, pc.p3_chain, 1, 1);
        merged_chains[1] = _get_merged_chain(nc.p3_chain, pc.p5_chain, 0, 1);
    }
    
    return merged_chains;
    
}

ChainOPs
MotifTreeMerger::_non_helix_merge(
    ChainEndPairMap const & nc,
    ChainEndPairMap const & pc) {
    
    ChainOPs merged_chains(2);
    ChainOP p3_chain = nc.p3_chain->subchain(0, nc.p3_chain->length()-1);
    ChainOP p5_chain = nc.p5_chain->subchain(1, nc.p5_chain->length());
    if     ( nc.is_hairpin() && pc.is_hairpin() ) {
        merged_chains[0] = _get_merged_chain(nc.p5_chain, pc.p3_chain, 1, 1);
    }
    else if( nc.is_hairpin() ) {
        auto p3_chain = pc.p3_chain->subchain(0, nc.p3_chain->length()-1);
        auto p5_chain = pc.p5_chain->subchain(1, nc.p5_chain->length());
        merged_chains[0] = _get_merged_hairpin(p3_chain, p5_chain, nc.p5_chain);
    }
    else if( pc.is_hairpin() ) {
        merged_chains[0] = _get_merged_hairpin(p5_chain, p3_chain, pc.p5_chain, 1);
        
    }
    else {
        merged_chains[0] = _get_merged_chain(p5_chain, pc.p3_chain, 1);
        merged_chains[1] = _get_merged_chain(p3_chain, pc.p5_chain);
    }
    
    return merged_chains;

}

ChainOP
MotifTreeMerger::_get_merged_chain(
    ChainOP const & c1,
    ChainOP const & c2,
    int join_by_3prime,
    int remove_overlap) {
    
    ResidueOPs chain1_res = c1->residues();
    ResidueOPs chain2_res = c2->residues();
    if(join_by_3prime) {
        std::reverse(chain1_res.begin(), chain1_res.end());
        std::reverse(chain2_res.begin(), chain2_res.end());
    }
    if(remove_overlap) { chain2_res.erase(chain2_res.begin()); }
    for(auto const & r : chain2_res) { chain1_res.push_back(r); }
    if(join_by_3prime) {
        std::reverse(chain1_res.begin(), chain1_res.end());
    }
    ChainOP merged_chain (new Chain(chain1_res));
    return merged_chain;
}

ChainOP
MotifTreeMerger::_get_merged_hairpin(
    ChainOP const & c1,
    ChainOP const & c2,
    ChainOP const & hairpin,
    int join_by_3prime,
    int remove_overlap) {
    
    ChainOP merged_chain = _get_merged_chain(c1, hairpin, join_by_3prime, remove_overlap);
    ChainOP merged_chain2 = _get_merged_chain(merged_chain, c2, join_by_3prime, remove_overlap);
    return merged_chain2;
}


