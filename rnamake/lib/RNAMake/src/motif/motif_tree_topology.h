//
//  motif_tree_topology.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_topology__
#define __RNAMake__motif_tree_topology__

#include <stdio.h>
#include "base/types.h"
#include "secondary_structure/ss_tree.h"
#include "data_structure/graph/graph.h"

enum class MotifTopologyType { BP_STEP, TWOWAY, NWAY, HAIRPIN, TCONTACT, START};

struct MotifTopologyBP {
    MotifTopologyBP(
        String nname,
        Strings nchains):
    name(nname),
    chains(nchains)
    {}
    
    String name;
    Strings chains;
};

struct MotifTopology {
    MotifTopology(
        MotifTopologyType const & ntype,
        Strings const & nseq,
        Strings const & nss,
        std::vector<MotifTopologyBP> nbps):
    type(ntype),
    seq(nseq),
    ss(nss),
    bps(nbps)
    {}

    MotifTopologyType type;
    Strings seq;
    Strings ss;
    std::vector<MotifTopologyBP> bps;
};

using MotifTopologyOP = std::shared_ptr<MotifTopology>;

class MotifTreeTopology {
public:
    MotifTreeTopology(SS_Tree const &);
    
private:
    VectorUP<SS_Node>
    get_bp_nodes(
        SS_Node const &,
        std::map<int, int> const & );
    
    MotifTopologyOP
    build_motif_topology_node(
        VectorUP<SS_Node> const &,
        SS_Tree const &,
        MotifTopologyType const &);
    
    
   
private:
    Graph<MotifTopologyOP> graph_;
    int bp_step_size_;
    
};
 

#endif /* defined(__RNAMake__motif_tree_topology__) */
















