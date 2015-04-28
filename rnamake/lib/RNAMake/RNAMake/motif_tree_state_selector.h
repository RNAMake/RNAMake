//
//  motif_tree_state_selector.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/21/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_selector__
#define __RNAMake__motif_tree_state_selector__

#include <stdio.h>
#include <map>
#include "types.h"
#include "motif_type.h"
#include "motif_tree_state_library.h"
#include "motif_tree_state_search_node.fwd.h"
#include "motif_tree_state_selector.fwd.h"

typedef std::pair<MotifTreeStateOP, int> MTSTypePair;
typedef std::vector<MTSTypePair> MTSTypePairs;
typedef std::shared_ptr<StringIntMap> StringIntMapOP;

class SelectorNode {
public:
    SelectorNode(
        MotifTreeStateLibraryOP const & mts_lib,
        int index,
        int max_uses = 1000000000,
        int required_uses = 0):
        mts_lib_ ( mts_lib ),
        index_ ( index ),
        connections_ ( SelectorNodeOPs() ),
        required_uses_(required_uses),
        max_uses_(max_uses) {}
    
    ~SelectorNode() {}
    
public:
    
    void
    add_connection(SelectorNodeOP const & node) {
        connections_.push_back(node);
    }
    
public:
    
    inline
    SelectorNodeOPs const &
    connections() { return connections_; }
    
    inline
    MotifTreeStateLibraryOP const &
    mts_lib() const { return mts_lib_; }
    
    inline
    int const & index() const { return index_; }
    
    inline
    int const & max_uses() const { return max_uses_; }
    
    inline
    int const & required_uses() const { return required_uses_; }
    
public:
    
    inline
    void max_uses(int nmax_uses) { max_uses_ = nmax_uses; }
    
    inline
    void required_uses(int nrequired_uses) { required_uses_ = nrequired_uses; }
    

private:
    MotifTreeStateLibraryOP mts_lib_;
    SelectorNodeOPs connections_;
    int index_, max_uses_, required_uses_;
};

class MotifTreeStateSelector {
public:
    MotifTreeStateSelector(
        MotifTreeStateLibraryOPs const &,
        String mode = "helix_flank");
    
    MotifTreeStateSelector(
        SelectorNodeOPs const & nodes):
    clash_lists_ ( std::map<String, StringIntMapOP>() ),
    mts_type_pairs_ ( MTSTypePairs() ) {
        nodes_ = nodes;
        _setup_clash_lists();
    }
    
    ~MotifTreeStateSelector() {}
    
public:
    MTSTypePairs const &
    get_children_mts(MotifTreeStateSearchNodeOP const &);
    
    int
    is_valid_solution(MotifTreeStateSearchNodeOP const &);
    
private:
    
    void
    _setup_clash_lists();
    
public:
    inline
    SelectorNodeOPs const &
    nodes() { return nodes_; }
    
private:
    SelectorNodeOPs nodes_;
    std::map<String, StringIntMapOP> clash_lists_;
    MTSTypePairs mts_type_pairs_;
    
};

MotifTreeStateSelector
default_selector(
    MotifTypes,
    String mode = "helix_flank");

typedef std::shared_ptr<MotifTreeStateSelector> MotifTreeStateSelectorOP;


#endif /* defined(__RNAMake__motif_tree_state_selector__) */
