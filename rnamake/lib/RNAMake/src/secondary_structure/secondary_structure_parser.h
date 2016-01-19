//
//  secondary_structure_parser.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_parser__
#define __RNAMake__secondary_structure_parser__

#include <cassert>
#include <stdio.h>
#include <map>

//RNAMake
#include "data_structure/graph/graph.h"
#include "secondary_structure/residue.h"
#include "secondary_structure/basepair.h"
#include "secondary_structure/structure.h"
#include "secondary_structure/motif.h"
#include "secondary_structure/pose.h"

namespace sstruct {
    

enum NodeType {
    UNPAIRED = 0,
    PAIRED   = 1
};
    
    
struct NodeData {
    NodeData() {}
    
    NodeData(
        ResidueOPs const & nresidues,
        NodeType const & ntype):
    residues(nresidues),
    type(ntype)
    {}
    
    ResidueOPs residues;
    NodeType type;
    
};
    
    
class SecondaryStructureChainGraph {
public:
    SecondaryStructureChainGraph():
    graph_(GraphStatic<NodeData>())
    {}
    
    ~SecondaryStructureChainGraph() {}
    
public:
    typedef typename GraphStatic<NodeData>::iterator iterator;
    typedef typename GraphStatic<NodeData>::const_iterator const_iterator;
    
    iterator begin() { return graph_.begin(); }
    iterator end()   { return graph_.end(); }
    
    const_iterator begin() const { return graph_.begin(); }
    const_iterator end()   const { return graph_.end(); }
    
public:
    size_t
    size() { return graph_.size(); }
    
    GraphNodeOPs<NodeData> const &
    nodes() { return graph_.nodes(); }

public:
    int
    add_chain(
        NodeData const & data,
        int parent_index = -1,
        int orphan = 0) {
        
        auto parent = graph_.last_node();
        if(parent_index != -1) {
            parent = graph_.get_node(parent_index);
        }
        if(parent == nullptr) {
            return graph_.add_data(data, -1, -1, -1, 3);
        }
        
        return graph_.add_data(data, parent_index, 1, 0, 3, orphan);
        
    }
    
    int
    get_node_by_res(
        ResidueOP const & res) {
        for(auto const & n : graph_) {
            for(auto const & r : n->data().residues) {
                if(r->uuid() == res->uuid()) { return n->index(); }
            }
        }
        return -1;
    }
  
    void
    pair_res(
        int n_i,
        int n_j) {
        
        assert(graph_.get_node(n_i)->data().type == NodeType::PAIRED &&
               "unpaired node is being paired");
        assert(graph_.get_node(n_j)->data().type == NodeType::PAIRED &&
               "unpaired node is being paired");
        graph_.connect(n_i, n_j, 2, 2);
    }
    
    
private:
    GraphStatic<NodeData> graph_;
    
};


typedef std::shared_ptr<GraphNode<NodeData>> SSNodeOP;
typedef std::shared_ptr<SecondaryStructureChainGraph> SecondaryStructureChainGraphOP;


class SecondaryStructureParser {
public:
    SecondaryStructureParser():
    seen_(std::map<SSNodeOP, int>())
    {}
    
    ~SecondaryStructureParser() {}

public:
    
    SecondaryStructureChainGraphOP
    parse(
        String const &,
        String const &);
    
    MotifOPs
    parse_to_motifs(
        String const &,
        String const &);
    
    MotifOP
    parse_to_motif(
        String const &,
        String const &);
    
    PoseOP
    parse_to_pose(
        String const &,
        String const &);

private:
    
    ResidueOP
    _previous_res(
        ResidueOP const & r) {
        
        int i = (int)(std::find(residues_.begin(), residues_.end(), r) - residues_.begin());
        if(i == 0) { return nullptr; }
        else       { return residues_[i-1]; }
        
    }
    
    int
    _start_of_chain(
        ResidueOP const & r) {
        
        for(auto const & c : structure_->chains()) {
            if(c->first() == r) { return 1; }
        }
        return 0;
    }
    
    ResidueOP
    _get_bracket_pair(
        ResidueOP const & r_start) {
        
        int bracket_count = 0;
        int start = 0;
        for(auto const & r : residues_) {
            if(r_start == r && !start) {
                start = 1;
            }
            else if(start) {}
            else { continue; }
            
            if(r->dot_bracket() == "(") { bracket_count += 1; }
            else if(r->dot_bracket() == ")") {
                bracket_count -= 1;
                if(bracket_count == 0) { return r; }
            }
            
        }
        
        throw SecondaryStructureException("cannot find pair in _get_bracket_pair");
        
    }
    
    MotifOP
    _generate_motif(
        SSNodeOP const &);
    
    SSNodeOP
    _walk_nodes(
        SSNodeOP const &);

    MotifOP
    _build_motif(
        StructureOP const &);
    
private:
    BasepairOPs pairs_;
    ResidueOPs residues_;
    StructureOP structure_;
    std::map<SSNodeOP, int> seen_;
    ChainOP chain_;
    
    
};

    
}
#endif /* defined(__RNAMake__secondary_structure_parser__) */
