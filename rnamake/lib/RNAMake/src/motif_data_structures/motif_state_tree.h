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
    cur_state(std::make_shared<MotifState>(nref_state->copy()))
    {}
    
    inline
    MSTNodeData(
        MSTNodeData const & ndata):
    ref_state(std::make_shared<MotifState>(*ndata.ref_state)),
    cur_state(std::make_shared<MotifState>(*ndata.cur_state))
    {}
    
public:
    inline
    BasepairStateOP const &
    get_end_state(
        String const & name) {
        return cur_state->get_end_state(name);
    }
    
    MotifStateOP ref_state, cur_state;
    
};

typedef TreeNodeOP<MSTNodeDataOP> MotifStateTreeNodeOP;

class MotifStateTree {
public:
    
    MotifStateTree() {
        tree_ = TreeStatic<MSTNodeDataOP>();
        aligner_ = MotifStateAligner();
        queue_ = std::queue<MotifStateTreeNodeOP>();
        options_ = Options("MotifTreeStateOptions");
        setup_options();
    }
    
    MotifStateTree(
        MotifTreeOP const & mt) {
        
        tree_ = TreeStatic<MSTNodeDataOP>();
        aligner_ = MotifStateAligner();
        queue_ = std::queue<MotifStateTreeNodeOP>();
        options_ = Options("MotifTreeStateOptions");
        setup_options();
        
        int i = -1;
        for(auto const & n : *mt) {
            i++;
            auto ms = RM::instance().motif_state(n->data()->name(),
                                                 n->data()->end_ids()[0],
                                                 n->data()->ends()[0]->name());
            
            if(i == 0) {
                add_state(ms);
            }
            
            else {
                int j = add_state(ms, n->parent_index(), n->parent_end_index());
                if(j == -1) {
                    throw MotifStateTreeException(
                        "could not convert motif tree to motif state tree");
                }
            }
            
        }
        
    }
    
    MotifStateTree(
        MotifStateTree const & mst):
    options_(Options("MotifStateTreeOptions")),
    tree_(TreeStatic<MSTNodeDataOP>(mst.tree_)) {
        
        for(auto const & n : mst) {
            tree_.get_node(n->index())->data() = std::make_shared<MSTNodeData>(*n->data());
        }
        
        options_ = Options(mst.options_);
        update_var_options();
    }
    
    
    ~MotifStateTree() {}
    
public: //iterators
    
    typedef typename TreeStatic<MSTNodeDataOP>::iterator iterator;
    typedef typename TreeStatic<MSTNodeDataOP>::const_iterator const_iterator;
    
    iterator begin() { return tree_.begin(); }
    iterator end()   { return tree_.end(); }
    
    const_iterator begin() const { return tree_.begin(); }
    const_iterator end()   const { return tree_.end(); }
    
public:
    
    int
    add_state(
        MotifStateOP const & state,
        int parent_index=-1,
        int parent_end_index=-1,
        String parent_end_name="",
        bool forced=false);
    
    int
    add_mst(
        MotifStateTreeOP const & mst,
        int parent_index=-1,
        int parent_end_index=-1,
        String parent_end_name="",
        bool forced=false);

    int
    setup_from_mt(
        MotifTreeOP const &);
    
    MotifTreeOP
    to_motif_tree();
    
    void
    replace_state(
        int i, 
        MotifStateOP const &);
    
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
    
    String
    topology_to_str();
    
public: //motif tree wrappers
    
    void
    write_pdbs(
        String const & fname = "nodes") {
        to_motif_tree()->write_pdbs(fname);
    }
    
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
    
    
private:
    void
    setup_options();
    
    void
    update_var_options();
    
    
private:
    TreeStatic<MSTNodeDataOP> tree_;
    std::queue<MotifStateTreeNodeOP> queue_;
    MotifStateAligner aligner_;
    MotifConnectionOPs connections_;
    Options options_;
    int sterics_;
    float clash_radius_;
    
};


#endif /* defined(__RNAMake__motif_state_tree__) */
