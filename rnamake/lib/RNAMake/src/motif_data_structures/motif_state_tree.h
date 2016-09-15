//
//  motif_state_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_tree__
#define __RNAMake__motif_state_tree__

#include <stdio.h>
#include <queue>

//RNAMake Headers
#include "base/types.h"
#include "base/option.h"
#include "data_structure/tree/tree.h"
#include "data_structure/tree/tree_node.h"
#include "motif/motif_state.h"
#include "motif/motif_state_aligner.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_state_tree.fwd.h"
#include "motif_data_structures/motif_connection.h"

class MotifStateTreeException : public std::runtime_error {
public:
    MotifStateTreeException(
        String const & message):
    std::runtime_error(message)
    {}
};

struct MSTNodeData {
public:
    inline
    MSTNodeData(
        MotifStateOP const & nref_state):
    ref_state(nref_state),
    cur_state(std::make_shared<MotifState>(*nref_state))
    {}
    
    inline
    MSTNodeData(
        MSTNodeData const & ndata):
    ref_state(std::make_shared<MotifState>(*ndata.ref_state)),
    cur_state(std::make_shared<MotifState>(*ndata.cur_state))
    {}
    
public: //wrappers for current state
    
    inline
    BasepairStateOP
    get_end_state(String const & name) { return cur_state->get_end_state(name); }
    
    inline
    int
    get_end_index(String const & name) { return cur_state->get_end_index(name); }
    
    inline
    String const &
    name() { return cur_state->name(); }
    
    inline
    int
    block_end_add() { return cur_state->block_end_add(); }
    
    inline
    String const &
    end_name(int i) { return cur_state->end_names()[i]; }
    
    inline
    Uuid const &
    uuid() { return cur_state->uuid(); }
    
public:
    MotifStateOP ref_state, cur_state;
    
};

typedef TreeNodeOP<MSTNodeDataOP> MotifStateTreeNodeOP;

class MotifStateTree {
public:
    
    MotifStateTree();
    
    MotifStateTree(MotifTreeOP const &);
    
    MotifStateTree(MotifStateTree const & );
    
    ~MotifStateTree() {}
    
public: //iterators
    
    typedef typename TreeStatic<MSTNodeDataOP>::iterator iterator;
    typedef typename TreeStatic<MSTNodeDataOP>::const_iterator const_iterator;
    
    iterator begin() { return tree_.begin(); }
    iterator end()   { return tree_.end(); }
    
    const_iterator begin() const { return tree_.begin(); }
    const_iterator end()   const { return tree_.end(); }
    
    
private: // add function helpers
    
    MotifStateTreeNodeOP
    _get_parent(int);
    
    Ints
    _get_available_parent_end_pos(
        MotifStateTreeNodeOP const &,
        int);
    
    int
    _get_parent_index_from_name(
        MotifStateTreeNodeOP const &,
            String const &);
    
    inline
    int
    _steric_clash(
        MSTNodeDataOP const & new_data) {
        
        float dist;
        for(auto const & n : tree_) {
            for(auto const & b1 : n->data()->cur_state->beads()) {
                for(auto const & b2 : new_data->cur_state->beads()) {
                    dist = b1.distance(b2);
                    if (dist < clash_radius_) { return 1; }
                }
            }
        }
        return 0;
    }
    
    int
    _get_connection_end(
        MotifStateTreeNodeOP const &,
        String const &);
    
public: // add functions
    
    int
    add_state(
        MotifStateOP const & state,
        int parent_index=-1,
        int parent_end_index=-1);
    
    int
    add_state(
        MotifStateOP const & state,
        int parent_index,
        String const & parent_end_name);
    
    int
    add_mst(
        MotifStateTreeOP const & mst,
        int parent_index=-1,
        int parent_end_index=-1);
    
    int
    add_mst(
        MotifStateTreeOP const &,
        int,
        String const &);
    
    
    void
    add_connection(
        int,
        int,
        String const &,
        String const &);
    
    void
    replace_state(
        int i,
        MotifStateOP const &);

public: // outputting functions
    
    MotifTreeOP
    to_motif_tree();
    
    String
    topology_to_str();
    
public: // getters
    
    inline
    Points
    centers() {
        auto centers = Points();
        for(auto const & n : tree_) {
            for(auto const & b : n->data()->cur_state->beads()) {
                centers.push_back(b);
            }
        }
        return centers;
    }
    
public: //motif tree wrappers
    
    void
    write_pdbs(
        String const & fname = "nodes") {
        to_motif_tree()->write_pdbs(fname);
    }
    
    inline
    RNAStructureOP 
    get_structure() { return to_motif_tree()->get_structure(); }
    
public: //tree wrapers
    size_t
    size() { return tree_.size(); }
    
    MotifStateTreeNodeOP const &
    last_node() { return tree_.last_node(); }
    
    MotifStateTreeNodeOP const &
    get_node(
        int i) {
        return tree_.get_node(i);
    }
    
    inline
    void
    increase_level() { tree_.increase_level(); }
    
    inline
    void
    decrease_level() { tree_.decrease_level(); }
    
    
public: //option wrappers
    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }
    
    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }
    
    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }
    
    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }
    
    
    template<typename T>
    void
    set_option_value(
        String const & name,
        T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }
    
private:
    void
    setup_options();
    
    void
    update_var_options();
    
    
private:
    TreeStatic<MSTNodeDataOP> tree_;
    std::queue<MotifStateTreeNodeOP> queue_;
    MotifStateAligner aligner_;
    MotifConnections connections_;
    Options options_;
    int sterics_;
    float clash_radius_;
    
};


#endif /* defined(__RNAMake__motif_state_tree__) */
