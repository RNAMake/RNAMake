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
#include "base/base.h"
#include "data_structure/graph/graph.h"
#include "data_structure/graph/graph_node.h"


#include "structure/residue_type_set.h"
#include "motif/motif.h"
#include "motif/motif_tree.fwd.h"
#include "motif/motif_tree_merger.h"
#include "motif/pose.h"



class MotifTree : public Base {
public:
    MotifTree();
    
    ~MotifTree() {}
    
public: //iterators
    
    typedef typename GraphStatic<MotifOP>::iterator iterator;
    typedef typename GraphStatic<MotifOP>::const_iterator const_iterator;
    
    iterator begin() { return graph_.begin(); }
    iterator end()   { return graph_.end(); }
    
    const_iterator begin() const { return graph_.begin(); }
    const_iterator end()   const { return graph_.end(); }
    
public:
    
    int
    add_motif(
        MotifOP const & m,
        int parent_index = -1,
        int parent_end_index = -1,
        String parent_end_name = "");
    
    void
    write_pdbs(
        String const & fname = "nodes");
    
    PoseOP
    to_pose();
    
    void
    to_pdb(
        String fname = "mt.pdb");
    
    void
    add_connection(
        int,
        int,
        String const &);
    
    void
    remove_node(
        int i = -1);
    
    void
    remove_node_level(int level=-1);
    
    String
    topology_to_str();
    

public: //graph wrappers
    
    inline
    void
    increase_level() {
        graph_.increase_level();
    }
    
    MotifTreeNodeOP
    get_node(int i) { return graph_.get_node(i); }
    
    size_t
    size() { return graph_.size(); }
    

    

public: //getters
    
    inline
    MotifTreeNodeOP const &
    last_node() { return graph_.last_node(); }
    
public: //setters
    
    
        
    
protected:
    
    void
    setup_options();
    
private:
    
    void
    update_var_options();
    
    int
    _steric_clash(
        MotifOP const &);
    

private:
    GraphStatic<MotifOP> graph_;
    MotifTreeMerger merger_;
    float clash_radius_;
    int level_;
    //options ... need a better way
    int sterics_;
    
    
};



MotifTree
str_to_motif_tree(
    String const &,
    ResidueTypeSet const &);



#endif /* defined(__RNAMake__motif_tree__) */
