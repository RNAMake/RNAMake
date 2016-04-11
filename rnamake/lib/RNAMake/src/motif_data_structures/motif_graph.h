//
//  motif_graph.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_graph__
#define __RNAMake__motif_graph__

#include <stdio.h>
#include <map>

//RNAMake Headers
#include "base/option.h"
#include "data_structure/graph/graph.h"
#include "motif/motif.h"
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_merger.h"

enum MotifGraphStringType {
    OLD,
    MG,
    TOP
};

class MotifGraphException : public std::runtime_error {
public:
    MotifGraphException(
        String const & message) :
    std::runtime_error(message)
    {}
};

class MotifGraph {
public:
    
    MotifGraph():
    graph_(GraphStatic<MotifOP>()),
    merger_(MotifMerger()),
    clash_radius_(2.5),
    sterics_(1),
    options_(Options("MotifGraphOptions")),
    aligned_(std::map<int, int>())
    {  setup_options(); }
    
    MotifGraph(String const &);
    
    MotifGraph(String const &,
               MotifGraphStringType const &);
    
    MotifGraph(
        MotifGraph const & mg):
    options_(Options("MotifGraphOptions")),
    graph_(GraphStatic<MotifOP>(mg.graph_)) {
        auto motifs = MotifOPs();
        // dear god this is horrible but cant figure out a better way to do a copy
        for(auto const & n : mg.graph_.nodes()) {
            graph_.get_node(n->index())->data() = std::make_shared<Motif>(*n->data());
            motifs.push_back(graph_.get_node(n->index())->data());
        }
        options_ = Options(mg.options_);
        merger_ = MotifMerger(mg.merger_, motifs);
        
    }

    ~MotifGraph() {}

public: //iterators
    
    typedef typename GraphStatic<MotifOP>::iterator iterator;
    typedef typename GraphStatic<MotifOP>::const_iterator const_iterator;
    
    iterator begin() { return graph_.begin(); }
    iterator end()   { return graph_.end(); }
    
    const_iterator begin() const { return graph_.begin(); }
    const_iterator end()   const { return graph_.end(); }
    
public:

    void
    remove_motif(int);
    
    void
    remove_level(int level);
    
    
    void
    replace_ideal_helices();
    
    void
    replace_helical_sequence(sstruct::PoseOP const &);
    
    
    sstruct::PoseOP
    designable_secondary_structure() {
        auto ss = merger_.secondary_structure();
        auto ss_r = sstruct::ResidueOP(nullptr);
        
        for(auto const & n : graph_) {
            if(n->data()->name() != "HELIX.IDEAL") { continue;}
            for(auto const & r : n->data()->residues()) {
                ss_r= ss->get_residue(r->uuid());
                if(ss_r != nullptr) {
                    ss_r->name("N");
                }
            }
        }
        
        return ss;
    }
    
    void
    write_pdbs(String const & fname = "nodes");
    
    inline
    Beads
    beads() {
        Beads beads;
        for(auto const & n : graph_.nodes()) {
            std::copy(n->data()->beads().begin(),
                      n->data()->beads().end(),
                      std::inserter(beads, beads.end()));
        }
        return beads;
    }
    
    String
    topology_to_str();
    
    String
    topology_to_str_new();

    
public: //add motif interface
    
    int
    add_motif(
        MotifOP const & m,
        int parent_index = -1,
        int parent_end_index = -1);
    
    int
    add_motif(
        String const & m_name,
        int parent_index = -1,
        int parent_end_index = 1);
    
    int
    add_motif(
        String const & m_name,
        String const & m_end_name,
        int parent_index = -1,
        int parent_end_index = 1);
    
    int
    add_motif(
        String const &,
        int,
        String const &);
    
    int
    add_motif(
        MotifOP const &,
        int,
        String const &);

public:
    
    void
    add_motif_tree(
        MotifTreeOP const & mt,
        int parent_index = -1);
    
    void
    add_motif_tree(
        MotifTreeOP const &,
        int,
        String const &);
    

    void
    add_connection(
        int,
        int,
        String const &);
    
    void
    add_connection(
        int,
        int,
        String const &,
        String const &);
    
    
    BasepairOP const &
    get_end(int);
    
    BasepairOP const &
    get_end(
        int,
        String const &);
    
    
    GraphNodeOPs<MotifOP>
    unaligned_nodes();
    
    inline
    GraphNodeOP<MotifOP> const &
    get_node_by_id(
        Uuid const & uuid) {
        for(auto const & n : graph_.nodes()) {
            if(n->data()->id() == uuid) {
                return n;
            }
        }
        throw std::runtime_error("could not find node by id");
    }

    
public: //Graph Wrappers
    inline
    size_t
    size() { return graph_.size(); }
    
    inline
    void
    increase_level() { return graph_.increase_level(); }
    
    inline
    GraphNodeOP<MotifOP> const &
    last_node() { return graph_.last_node(); }
    
    inline
    GraphNodeOP<MotifOP>
    oldest_node() { return graph_.oldest_node(); }
    
    inline
    GraphNodeOP<MotifOP> const &
    get_node(int i) { return graph_.get_node(i); }
    


public: //Motif Merger Wrappers
    
    inline
    RNAStructureOP const &
    get_structure() {
        return merger_.get_structure();
    }
    
    
    inline
    void
    to_pdb(
        String const fname = "test.pdb",
        int renumber = -1) {
        return merger_.to_pdb(fname, renumber);
    }
    
    sstruct::PoseOP
    secondary_structure() {
        return merger_.secondary_structure();
    }
    
    
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
    int
    _steric_clash(
        MotifOP const &);
    
    int
    _add_motif_to_graph(
        MotifOP const & m,
        int parent_index = -1,
        int parent_end_index = -1);
    
private:
    void
    setup_options();
    
    void
    update_var_options();
    
private:
    void
    _setup_from_top_str(String const &);
    
    void
    _setup_from_str(String const &);
    
private:
    GraphStatic<MotifOP> graph_;
    MotifMerger merger_;
    Options options_;
    std::map<int, int> aligned_;
    //options
    float clash_radius_;
    bool sterics_;

};

typedef std::shared_ptr<MotifGraph> MotifGraphOP;

#endif /* defined(__RNAMake__motif_graph__) */
