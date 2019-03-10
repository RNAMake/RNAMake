//
//  motif_state_graph.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef motif_state_graph_hpp
#define motif_state_graph_hpp

#include <stdio.h>
#include "base/option.h"
#include "data_structure/graph/graph.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_data_structures/motif_state_node.hpp"

class MotifStateGraphException : public std::runtime_error {
public:
    MotifStateGraphException(
        String const & message):
    std::runtime_error(message)
    {}
};

class MotifStateGraph {
public:
    MotifStateGraph();
    
    MotifStateGraph(MotifGraphOP const &);
    
    MotifStateGraph(MotifStateGraph const &);
    
    ~MotifStateGraph() {}

public: //iterators
    
    typedef typename GraphStatic<MSNodeDataOP>::iterator iterator;
    typedef typename GraphStatic<MSNodeDataOP>::const_iterator const_iterator;
    
    iterator begin() { return graph_.begin(); }
    iterator end()   { return graph_.end(); }
    
    const_iterator begin() const { return graph_.begin(); }
    const_iterator end()   const { return graph_.end(); }

public:

    typedef typename GraphNodeOPs<MSNodeDataOP>::iterator node_iterator;
    typedef typename GraphNodeOPs<MSNodeDataOP>::const_iterator node_const_iterator;

    node_iterator node_begin() { _update_align_list(); return align_list_.begin(); }
    node_iterator node_end()   { return align_list_.end(); }



private:
    void
    _setup_from_mg(MotifGraphOP const &);
    
public:
    size_t
    size() { return graph_.size(); }
    
    
private://add function helpers
    
    GraphNodeOP<MSNodeDataOP>
    _get_parent(
        String const &,
        int);

    Ints
    _get_available_parent_end_pos(
        GraphNodeOP<MSNodeDataOP> const &,
        int);
    
    int
    _get_parent_index_from_name(
        GraphNodeOP<MSNodeDataOP> const &,
        String const &);

    
    int
    _get_connection_end(
        GraphNodeOP<MSNodeDataOP> const &,
        String const &);
    
    inline
    int
    _steric_clash(
        MSNodeDataOP const & new_data) {
        float dist;
        for(auto const & n : graph_.nodes()) {
            for(auto const & b1 : n->data()->cur_state->beads()) {
                for(auto const & b2 : new_data->cur_state->beads()) {
                    dist = b1.distance(b2);
                    if (dist < clash_radius_) { return 1; }
                }
            }
        }
        return 0;
    }
    
public: // add functions
    
    int
    add_state(
            MotifStateOP const & state,
            int parent_index=-1,
            int parent_end_index=-1,
            int orphan = 0);
    
    int
    add_state(
        MotifStateOP const & state,
        int parent_index,
        String const & parent_end_name);
    
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

    void
    remove_state(int);

public: //remove functions

    void
    remove_level(int level);
    
    
public: // graph wrappers

    inline
    GraphNodeOP<MSNodeDataOP>
    last_node() { return graph_.last_node(); }

    inline
    GraphNodeOP<MSNodeDataOP> const &
    get_node(int i) const { return graph_.get_node(i); }
    
    inline
    GraphNodeOP<MSNodeDataOP> const
    get_node(Uuid const & uuid) const {
        for(auto const & n : graph_) {
            if(n->data()->uuid() == uuid) {
                return n;
            }
        }
        throw MotifStateGraphException("cannot get node with uuid no motif has it in this tree");
    }
    
    inline
    GraphNodeOP<MSNodeDataOP> const
    get_node(String const & m_name) const {
        auto node = GraphNodeOP<MSNodeDataOP>(nullptr);
        for(auto const & n : graph_) {
            if(n->data()->name() == m_name) {
                if(node != nullptr) {
                    throw MotifStateGraphException(
                            "cannot get node with name: " + m_name + " there is more then one motif "
                            "with this name");
                }
                
                node = n;
            }
        }
        
        if(node == nullptr) {
            throw MotifStateGraphException(
                    "cannot get node with name: " + m_name + " there is no motif in the tree with "
                    "this name");
        }
        
        return node;
    }

    void
    increase_level() { return graph_.increase_level(); }
    
    
public: // motif graph wrappers
    
    MotifGraphOP
    to_motif_graph();
    
    
private: // misc functions
    
    void
    _update_align_list();
    
    void
    _align_states(int pos=-1);

    
public: // getters
    
    GraphNodeOPs<MSNodeDataOP> const
    unaligned_nodes() const;

    
public: //Options Wrappers
    
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
    GraphStatic<MSNodeDataOP> graph_;
    base::Options options_;
    std::map<int, int> aligned_;
    GraphNodeOPs<MSNodeDataOP> align_list_;
    int update_align_list_;
    //options
    float clash_radius_;
    bool sterics_;
};

typedef std::shared_ptr<MotifStateGraph> MotifStateGraphOP;

#endif /* motif_state_graph_hpp */
