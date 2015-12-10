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

//RNAMake Headers
#include "data_structure/graph/graph.h"
#include "motif/motif.h"
#include "motif_data_structures/motif_merger.h"

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
    sterics_(1)
    {}
    
    MotifGraph(
        MotifGraph const & mg):
    graph_(GraphStatic<MotifOP>(mg.graph_)) {
        auto motifs = MotifOPs();
        // dear god this is horrible but cant figure out a better way to do a copy
        for(auto const & n : mg.graph_.nodes()) {
            graph_.get_node(n->index())->data() = std::make_shared<Motif>(*n->data());
            motifs.push_back(graph_.get_node(n->index())->data());
        }
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

    size_t
    size() { return graph_.size(); }
    
    void
    remove_motif(int);
    
    void
    replace_ideal_helices();
    
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

public:
    
    inline
    RNAStructureOP const &
    get_structure() {
        return merger_.get_structure();
    }
    
    sstruct::PoseOP
    secondary_structure() {
        return merger_.secondary_structure();
    }
    
    void
    write_pdbs(String const & fname = "nodes");
    

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
    GraphStatic<MotifOP> graph_;
    MotifMerger merger_;
    float clash_radius_;
    int sterics_;

};

typedef std::shared_ptr<MotifGraph> MotifGraphOP;

#endif /* defined(__RNAMake__motif_graph__) */
