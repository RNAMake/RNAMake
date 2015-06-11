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
        Strings nchains,
        Ints nid):
    name(nname),
    chains(nchains),
    id(nid) {}
    
    bool
    operator == (MotifTopologyBP const & rhs) const {
        if(id[0] == rhs.id[0] && id[1] == rhs.id[1]) {  return true; }
        if(id[1] == rhs.id[0] && id[0] == rhs.id[1]) {  return true; }
        return false;
    }
    
    String name;
    Strings chains;
    Ints id;
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
    ends(nbps) {
        name = "";
        int i = 0;
        for(int i = 0; i < seq.size(); i++) {
            name += seq[i] + "_";
            for(auto const & e : ss[i]) {
                if     (e == '(') { name += "L"; }
                else if(e == ')') { name += "R"; }
                else              { name += "U"; }
            }
            if(i != (int)seq.size()-1) { name += "_"; }
        }
    
    }

    MotifTopologyType type;
    Strings seq;
    Strings ss;
    std::vector<MotifTopologyBP> ends;
    String name;
};

using MotifTopologyOP = std::shared_ptr<MotifTopology>;
using MTT_Node = GraphNodeOP<MotifTopologyOP>;
using MTT_Nodes = std::vector<MTT_Node>;

struct ConnectivityInfo {
    ConnectivityInfo() {}
    
    int n_index;
    int i;
    int j;
};

class MotifTreeTopology {
public:
    MotifTreeTopology(SS_Tree const &);
    
public:
    
    inline
    auto nodes() {
        return graph_.nodes();
    }
    
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
    
    ConnectivityInfo
    get_connectivity_info(
        MotifTopologyOP const & n);
    
private:
    GraphStatic<MotifTopologyOP> graph_;
    int bp_step_size_;
    
};
 

#endif /* defined(__RNAMake__motif_tree_topology__) */
















