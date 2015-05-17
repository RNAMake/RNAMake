//
//  motif_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree__
#define __RNAMake__motif_tree__

#include <stdio.h>
#include <map>

//RNAMake Headers
#include "base/types.h"
#include "structure/residue_type_set.h"
#include "motif/motif_tree.fwd.h"
#include "motif/motif_tree_node.h"
#include "motif/motif.h"
#include "motif/motif_tree_merger.h"
#include "motif/pose.h"


class MotifTree {
public:
    MotifTree();
    MotifTree(
        MotifOP const &);
    
    MotifTree(
        String const &,
        ResidueTypeSet const &);
    
    ~MotifTree() {
        last_node_ = nullptr;
        for(auto const & n : nodes_) {
            for(auto & c : n->connections()) { c->disconnect(); }
        }
     }
    
public:
    
    MotifTreeNodeOP
    add_motif(
        MotifOP const & m,
        MotifTreeNodeOP parent = nullptr,
        int end_index = -1,
        int parent_end = -1,
        int end_flip = -1);
    
    void
    write_pdbs(
        String const & fname = "nodes");
    
    PoseOP
    to_pose(
        int include_head = 0);
    
    void
    to_pdb(
        String fname = "mt.pdb",
        int include_head = 0);
    
    void
    _add_connection(
        MotifTreeNodeOP const &,
        MotifTreeNodeOP const &,
        float cutoff = 25);
    
    inline
    void
    increase_level() { level_ += 1; }
    
    void
    remove_node(
        MotifTreeNodeOP const &);
    
    void
    inline
    remove_last_node() { remove_node(last_node_); }
    
    void
    remove_node_level(int level=-1);

public: //getters
    
    inline
    MotifTreeNodeOPs const &
    nodes() const { return nodes_; }
    
    inline
    MotifTreeNodeOP const &
    last_node() { return last_node_; }
    
public: //setters
    
    inline
    void
    sterics(int nsterics) { sterics_ = nsterics; }
    
private:
    
    MotifTreeNodeOP
    _add_motif(
        MotifTreeNodeOP const &,
        MotifOP const &,
        BasepairOP const &,
        BasepairOP const &,
        Ints const &);
    
    int
    _steric_clash(
        MotifOP const &,
        BasepairOP const &);
    
    void
    _update_beads(
        MotifTreeNodeOP const &,
        MotifTreeNodeOP const &);

private:
    MotifTreeNodeOPs nodes_;
    MotifTreeNodeOP last_node_;
    MotifTreeMerger merger_;
    float clash_radius_;
    int level_;
    //options ... need a better way
    int sterics_, full_beads_first_res_;
    
    
};



MotifTree
str_to_motif_tree(
    String const &,
    ResidueTypeSet const &);



#endif /* defined(__RNAMake__motif_tree__) */
