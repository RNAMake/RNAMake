//
//  motif_toplogy.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_topology__
#define __RNAMake__motif_topology__

#include <stdio.h>

#include "motif_data_structure/motif_tree.h"
#include "motif_data_structure/motif_graph.h"

namespace motif_data_structure {

class MotifTopologyException : public std::runtime_error {
public:
    MotifTopologyException(String const & message) :
            std::runtime_error(message) {}
};


class GraphtoTree {
public:
    GraphtoTree() {}

    ~GraphtoTree() {}

private:
    struct _GraphtoTreeNode {
        inline
        _GraphtoTreeNode(
                data_structure::graph::GraphNodeOP<motif::MotifOP> const & nparent,
                int nparent_end_index,
                data_structure::graph::GraphNodeOP<motif::MotifOP> const & nnode,
                motif::MotifOP nmotif) :
                parent(nparent),
                parent_end_index(nparent_end_index),
                node(nnode),
                motif(nmotif) {}

        data_structure::graph::GraphNodeOP<motif::MotifOP> parent;
        int parent_end_index;
        data_structure::graph::GraphNodeOP<motif::MotifOP> node;
        motif::MotifOP motif;
    };

    typedef std::shared_ptr<_GraphtoTreeNode> _GraphtoTreeNodeOP;
    typedef std::vector<_GraphtoTreeNodeOP> _GraphtoTreeNodeOPs;


public:
    MotifTreeOP
    convert(
            MotifGraphOP const & mg,
            data_structure::graph::GraphNodeOP<motif::MotifOP> start = nullptr,
            int start_end_index = -1,
            data_structure::graph::GraphNodeOP<motif::MotifOP> last_node = nullptr);

private:
    _GraphtoTreeNodeOP
    _get_start_node(
            MotifGraphOP const &,
            data_structure::graph::GraphNodeOP<motif::MotifOP> const &,
            int);

    motif::MotifOP
    _get_reoriented_motif(
            motif::MotifOP const &,
            int);

    _GraphtoTreeNodeOPs
    _get_new_nodes(
            _GraphtoTreeNodeOP const &);

    int
    _get_new_parent_end_index(
            data_structure::graph::GraphNodeOP<motif::MotifOP> const &,
            data_structure::graph::GraphConnectionOP<motif::MotifOP> const &);

private:
    MotifTreeOP mt_;
};

MotifTreeOP
graph_to_tree(
        MotifGraphOP const & mg,
        data_structure::graph::GraphNodeOP<motif::MotifOP> start = nullptr,
        structure::BasepairOP last_end = nullptr);

}


#endif /* defined(__RNAMake__motif_toplogy__) */
