//
//  motif_tree_merger.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_merger.h"
#include "motif_tree.h"
#include "chain.h"
#include "structure.h"

PoseOP
MotifTreeMerger::merge(MotifTree const & mt, int include_head) {
    include_head_ = include_head;
    if(mt.nodes().size() == 2 && include_head_ == 0) {
        MotifOP m (new Motif( mt.nodes()[1]->motif()->copy() ));
        return PoseOP ( new Pose(m) );
    }

    for(auto const & n : mt.nodes())  { nodes_.push_back(n); }
    
    int i = -1;
    for( auto const & n : nodes_) {
        i++;
        if(i == 0 && include_head_ == 0) { continue; }
        for( auto const & c : n->motif()->chains()) {
            chains_.push_back(c->subchain(0, (int)c->residues().size()));
        }
    }
    
    MotifTreeNodeOP start_node = nodes_[0];
    if(include_head_ == 0) { start_node = nodes_[1]; }
    _merge_chains_in_node(start_node);
    PoseOP new_pose = _build_pose();
    
    seen_connections_.clear(); //= std::map<MotifTreeConnectionOP, int>();
    chains_.resize(0);
    nodes_.resize(0);
    
    return new_pose;
}

PoseOP
MotifTreeMerger::_build_pose() {
    StructureOP new_structure ( new Structure() );
    ChainOPs new_chains;
        
    for( auto const & c : chains_) {
        new_chains.push_back( ChainOP( new Chain(c->copy())) );
    }
    new_structure->chains(new_chains);
    new_structure->renumber();
    new_structure->_cache_coords();
    ResidueOPs residues = new_structure->residues();
    std::map<String, ResidueOP> uuids;
    for(auto const & res : residues) { uuids[res->uuid().s_uuid() ] = res; }
    BasepairOPs basepairs;
    std::map<String, int> designable;
    for(auto const & node : nodes_) {
        for(auto const & bp : node->motif()->basepairs()) {
            if(uuids.find(bp->res1()->uuid().s_uuid()) ==  uuids.end()) { continue; }
            if(uuids.find(bp->res2()->uuid().s_uuid()) ==  uuids.end()) { continue; }
            BasepairOP cbp (new Basepair(bp->copy()));
            cbp->res1( uuids[bp->res1()->uuid().s_uuid()]);
            cbp->res2( uuids[bp->res2()->uuid().s_uuid()]);
            if(node->motif()->mtype() == HELIX) { designable[cbp->uuid().s_uuid()] = 1; }
            basepairs.push_back(cbp);
        }
    }

    PoseOP new_pose (new Pose(new_structure, basepairs));
    new_pose->designable(designable);
    new_pose->_cache_basepair_frames();
    
    
    return new_pose;
    
}

void
MotifTreeMerger::_merge_chains_in_node(
    MotifTreeNodeOP const & node) {
    for(auto const & c : node->connections() ) {
        if(seen_connections_.find(c) != seen_connections_.end() ) { continue; }
        seen_connections_[c] = 1;
        MotifTreeNodeOP partner = c->partner(node);
        if(include_head_ == 0 && partner == nodes_[0]) { continue; }
        BasepairOP end  = c->motif_end(node);
        BasepairOP pend = c->motif_end(partner);
        ChainEndPairMap node_chains = _find_chains_for_end(end);
        ChainEndPairMap partner_chains = _find_chains_for_end(pend);
        ChainOPs merged_chains;
        if(partner->motif()->mtype() == HELIX) { merged_chains = _helix_merge(node_chains, partner_chains); }
        else                                   { merged_chains = _non_helix_merge(node_chains, partner_chains);}
        
        std::map<ChainOP, int> used_chains;
        ChainOPs new_chains;
        for (auto const & chain : node_chains.chains()) { used_chains[chain] = 1; }
        for (auto const & chain : partner_chains.chains()) { used_chains[chain] = 1; }
        for (auto const & chain : chains_) {
            if (used_chains.find(chain) == used_chains.end()) { new_chains.push_back(chain); }
        }
        for (auto const & chain : merged_chains) {
            if (chain != NULL) { new_chains.push_back(chain); }
        }
        chains_ = new_chains;
        _merge_chains_in_node(partner);
    }
    
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
        throw "Could not find chain for end";
    }
    
    if(chain_infos[0].pos == 0) { return ChainEndPairMap(chain_infos[0].chain, chain_infos[1].chain); }
    else {                        return ChainEndPairMap(chain_infos[1].chain, chain_infos[0].chain); }
    
}

ChainOPs
MotifTreeMerger::_helix_merge(
    ChainEndPairMap const & nc,
    ChainEndPairMap const & pc) {
    
    ChainOPs merged_chains(2);
    if     ( nc.is_hairpin() && pc.is_hairpin() ) {
        
        merged_chains[0] = _get_merged_chain(nc.p5_chain, pc.p3_chain, 1, 1);

        //throw "cannot merge an hairpin with another hairpin";
    
    }
    else if( nc.is_hairpin() ) {
        ChainOP p3_chain = pc.p3_chain->subchain(0, (int)pc.p5_chain->residues().size()-1);
        ChainOP p5_chain = pc.p5_chain->subchain(1, (int)pc.p5_chain->residues().size());
        merged_chains[0] = _get_merged_hairpin(p3_chain, p5_chain, nc.p5_chain);
    }
    else if( pc.is_hairpin() ) {
        merged_chains[0] = _get_merged_hairpin(nc.p5_chain, nc.p3_chain, pc.p5_chain, 1, 1);
        
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
    ChainOP p3_chain = nc.p3_chain->subchain(0, (int)nc.p3_chain->residues().size()-1);
    ChainOP p5_chain = nc.p5_chain->subchain(1, (int)nc.p5_chain->residues().size());
    if     ( nc.is_hairpin() && pc.is_hairpin() ) {
        merged_chains[0] = _get_merged_chain(nc.p5_chain, pc.p3_chain, 1, 1);
    }
    else if( nc.is_hairpin() ) {
        merged_chains[0] = _get_merged_hairpin(pc.p3_chain, pc.p5_chain, p5_chain);
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
    
    ChainOP merged_chain (new Chain());
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
    merged_chain->residues(chain1_res);
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


