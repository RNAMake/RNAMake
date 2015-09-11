//
//  ss_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_tree__
#define __RNAMake__ss_tree__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "data_structure/tree/tree.h"
#include "secondary_structure/secondary_structure.h"
#include "secondary_structure/ss_tree_node.h"
#include "secondary_structure/ss_tree_node.fwd.h"

namespace sstruct {

typedef TreeNodeOP<SS_NodeDataOP>  SS_TreeNodeOP;
typedef std::vector<SS_TreeNodeOP> SS_TreeNodeOPs;

//using SS_Nodes = std::vector<SS_Node>;

class SS_Tree {
public:
    SS_Tree(
        String const &,
        String const &);
    
    ~SS_Tree() {}
    
public:
    
    typedef typename TreeDynamic<SS_NodeDataOP>::iterator iterator;
    typedef typename TreeDynamic<SS_NodeDataOP>::const_iterator const_iterator;
    
    iterator begin() { return tree_.begin(); }
    iterator end()   { return tree_.end(); }
    
    const_iterator begin() const { return tree_.begin(); }
    const_iterator end()   const { return tree_.end(); }
    

public:
    inline
    int
    size() { return (int)tree_.size(); }
    
private:
    void
    _build_tree();
    
    std::vector<SS_NodeDataOP>
    _build_tree_level(
        int,
        int);
  
    SS_NodeDataOP
    _assign_new_node(
        int,
        int);
    
    int
    _get_brack_pair(
        int);
    
    int
    _get_dot_bounds(
        int,
        int);
    
    int
    _map_back_to_index(
        ResidueOP const &);
    
    SS_NodeDataOP
    _check_from_chain_ends(
        int,
        int);

    int
    _is_res_end_of_chain(
        int);
    
public: //getters
    SecondaryStructure const &
    secondary_structure() { return ss_; }
    
private:
    TreeDynamic<SS_NodeDataOP> tree_;
    ResidueOPs residues_; 
    SecondaryStructure ss_;
    
};
    
} //sstruct

#endif /* defined(__RNAMake__ss_tree__) */
