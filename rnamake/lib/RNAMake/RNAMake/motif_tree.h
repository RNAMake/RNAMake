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
#include "types.h"
#include "motif_tree.fwd.h"
#include "motif_tree_node.h"
#include "motif.h"
#include "uuid.h"
#include "residue_type_set.h"
#include "motif_tree_merger.h"
#include "pose.h"


class MotifTree {
public:
    MotifTree();
    MotifTree(MotifOP const &);
    MotifTree(String const &, ResidueTypeSet const &);
    ~MotifTree() {}
    
public:
    
    MotifTreeNodeOP
    add_motif(
        MotifOP const & m,
        MotifTreeNodeOP parent = NULL,
        int end_index = -1,
        int parent_end = -1,
        int end_flip = -1);
    
    void
    write_pdbs(String const & fname = "nodes");
    
    PoseOP
    to_pose();

public: //getters
    
    inline
    MotifTreeNodeOPs const &
    nodes() const { return nodes_; }
    
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
    //options ... need a better way
    int sterics_, full_beads_first_res_;
    
    
};



MotifTree
str_to_motif_tree(
    String const &,
    ResidueTypeSet const &);



#endif /* defined(__RNAMake__motif_tree__) */
