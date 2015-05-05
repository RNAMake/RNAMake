//
//  motif_tree_merger.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_merger__
#define __RNAMake__motif_tree_merger__

#include <stdio.h>
#include <map>

//RNAMake Headers
#include "structure/chain.fwd.h"
#include "motif/motif_tree.fwd.h"
#include "motif/motif_tree_node.h"
#include "motif/pose.h"

struct ChainEndPairMap {
public:
    ChainEndPairMap(
        ChainOP const & np5_chain,
        ChainOP const & np3_chain):
        p5_chain ( np5_chain),
        p3_chain ( np3_chain) {}
    
    ~ChainEndPairMap() {}

public:
    
    int
    is_hairpin() const {
        if(p5_chain == p3_chain) { return 1; }
        else                     { return 0; }
    }
    
    ChainOPs
    chains() const {
        ChainOPs chains(2);
        chains[0] = p5_chain;
        chains[1] = p3_chain;
        return chains;
    }
    
public:
    ChainOP p5_chain, p3_chain;
};

struct ChainInfo {
public:
    ChainInfo(
        ChainOP const & nchain,
        int npos,
        int nindex):
        chain ( nchain ),
        pos ( npos ),
        index ( nindex ) {}
    
public:
    ChainOP chain;
    int pos, index;
};

class MotifTreeMerger {
public:
    MotifTreeMerger() {
        seen_connections_ = std::map<MotifTreeConnectionOP, int>();
        chains_ = ChainOPs();
        nodes_ = MotifTreeNodeOPs();
        include_head_ = 0;
    }
    
public:
    PoseOP
    merge(MotifTree const &, int include_head = 0);
    
private:
    
    void
    _merge_chains_in_node(
        MotifTreeNodeOP const &);
    
    ChainEndPairMap
    _find_chains_for_end(
        BasepairOP const &);
    
    ChainOPs
    _helix_merge(
        ChainEndPairMap const &,
        ChainEndPairMap const &);
    
    ChainOPs
    _non_helix_merge(
        ChainEndPairMap const &,
        ChainEndPairMap const &);
    
    ChainOP
    _get_merged_chain(
        ChainOP const &,
        ChainOP const &,
        int join_by_3prime = 0,
        int remove_overlap = 0);
    
    ChainOP
    _get_merged_hairpin(
        ChainOP const &,
        ChainOP const &,
        ChainOP const &,
        int join_by_3prime = 0,
        int remove_overlap = 0);
    
    PoseOP
    _build_pose();
    
private:
    std::map<MotifTreeConnectionOP, int> seen_connections_;
    ChainOPs chains_;
    MotifTreeNodeOPs nodes_;
    int include_head_;
};

#endif /* defined(__RNAMake__motif_tree_merger__) */
