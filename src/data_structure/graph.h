//
// Created by Joseph Yesselman on 1/14/18.
//

#ifndef RNAMAKE_NEW_NEW_GRAPH_H
#define RNAMAKE_NEW_NEW_GRAPH_H

#include <queue>

#include <base/types.hpp>
#include <base/assertions.h>
#include <base/vector_container.h>
#include <data_structure/graph_base.h>
#include <data_structure/graph_adjacency_list.h>
#include <data_structure/graph_iter_list.h>

namespace data_structure {

template<typename DataType, typename AdjacencyListType, typename IterListType>
class _Graph {
public:
    _Graph() :
            adjacency_list_(AdjacencyListType()),
            iter_list_(IterListType()){}

    virtual
    ~_Graph() {}

public:
    typedef typename IterListType::const_iterator const_iterator;
    typedef typename IterListType::iterator iterator;

    iterator begin(){ return iter_list_.begin(); }
    iterator end()  { return iter_list_.end(); }

    const_iterator begin() const noexcept { return iter_list_.begin(); }
    const_iterator end() const noexcept   { return iter_list_.end(); }

public:
    void
    setup_transversal(
            Index start_n) {
        iter_list_.transversal(adjacency_list_, start_n);
    }

    void
    setup_path_transversal(
            Index start_n,
            Index end_n) {
        iter_list_.path_transversal(adjacency_list_, start_n, end_n);
    }

public:
    virtual
    inline
    Index
    add_node(
            DataType const & d,
            Size n_edges) {
        return adjacency_list_.add_node(d, n_edges);
    }

    inline
    void
    add_edge(
            NodeIndexandEdge const & nie1,
            NodeIndexandEdge const & nie2) {
        return adjacency_list_.add_edge(nie1, nie2);
    }

    inline
    void
    remove_node(
            Index ni) {
        return adjacency_list_.remove_node(ni);
    }

    inline
    void
    remove_edge(
            NodeIndexandEdge const & nie1,
            NodeIndexandEdge const & nie2) {
        return adjacency_list_.remove_edge(nie1, nie2);
    }

public:
    inline
    size_t
    get_num_nodes() const {
        return adjacency_list_.get_num_nodes();
    }

    inline
    size_t
    get_num_edges() const {
        return adjacency_list_.get_num_edges();
    }

    inline
    std::vector<Edge const *> const &
    get_node_edges(
            Index ni) const {
        return adjacency_list_.get_node_edges(ni);
    }

    inline
    Node<DataType> const &
    get_node(
            Index ni) const {
        return adjacency_list_.get_node(ni);
    }

    inline
    DataType const &
    get_node_data(
            Index ni) const {
        return adjacency_list_.get_node_data(ni);
    }

    inline
    DataType &
    get_node_data(
            Index ni) {
        return adjacency_list_.get_node_data(ni);
    }

    inline
    NodeIndexandEdge
    get_connected_node_info(
            NodeIndexandEdge const & nei) const {
        return adjacency_list_.get_connected_node_info(nei);
    }

    inline
    bool
    edge_between_nodes(
            Index n1,
            Index n2) const {
        return adjacency_list_.edge_between_nodes(n1, n2);
    }

    inline
    bool
    edge_index_empty(
            Index ni,
            Index ei) const {
        return adjacency_list_.edge_index_empty(ni, ei);
    }


protected:
    AdjacencyListType adjacency_list_;
    mutable IterListType iter_list_; // needs to update to iterate

};

template <typename Datatype, typename EdgeType>
using _UndirectedGraph = _Graph<
        Datatype,
        AdjacencyList<Datatype, EdgeType>,
        IterList<Datatype, AdjacencyList<Datatype, EdgeType>>>;

template<typename DataType, typename EdgeType>
using _DirectedGraph = _Graph<
        DataType,
        DirectedAdjacencyList<DataType, EdgeType>,
        DirectedIterList<DataType, DirectedAdjacencyList<DataType, EdgeType>>>;

template<typename DataType, typename EdgeType>
class UndirectedGraph : public _UndirectedGraph <DataType, EdgeType> {
public:
    typedef _UndirectedGraph<DataType, EdgeType> BaseClass;
public:
    UndirectedGraph():
            BaseClass() {}

    UndirectedGraph(
            UndirectedGraph const & g) {
        this->adjacency_list_ = g.adjacency_list_;

    }

};

template<typename DataType>
using FixedEdgeUndirectedGraph = UndirectedGraph<DataType, FixedEdges>;

template<typename DataType>
using DynamicEdgedUndirectedGraph = UndirectedGraph<DataType, DynamicEdges>;


template<typename DataType, typename EdgeType>
class DirectedGraph : public _DirectedGraph<DataType, EdgeType> {
public:
    typedef _DirectedGraph<DataType, EdgeType> BaseClass;

public:
    DirectedGraph() :
            BaseClass() {}

    DirectedGraph(
            DirectedGraph const & g):
            BaseClass() {
        this->adjacency_list_ = g.adjacency_list_;
    }

public:
    inline
    void
    setup_sub_graph_transversal(
            Index start_n,
            Index end_n) {
        return this->iter_list_.sub_graph_transversal(this->adjacency_list_, start_n, end_n);
    }

public:
    inline
    Index
    add_node(
            DataType const & d,
            Size n_edges,
            Index n_end_index,
            NodeIndexandEdge const & pie) {
        return this->adjacency_list_.add_node(d, n_edges, n_end_index, pie);
    }

    inline
    Index
    add_node(
            DataType const & d,
            Size n_edges) {
        return this->adjacency_list_.add_node(d, n_edges);
    }

public:
    inline
    bool
    has_parent(
            Index ni) const {
        return this->adjacency_list_.has_parent(ni);
    }

    inline
    Index
    get_parent_index(
            Index ni) const {
        return this->adjacency_list_.get_parent_index(ni);
    }

    inline
    Index
    get_parent_end_index(
            Index ni) const {
        return this->adjacency_list_.get_parent_end_index(ni);
    }

public:
    Indexes
    get_root_indexes() {
        auto roots = Indexes();
        for(auto const & kv : this->adjacency_list_) {
            if(!this->adjacency_list_.has_parent(kv.first)) {
                roots.push_back(kv.first);
            }
        }
        return roots;
    }
};

template<typename DataType>
using FixedEdgeDirectedGraph = DirectedGraph<DataType, FixedEdges>;

template<typename DataType>
using DynamicEdgedDirectedGraph = DirectedGraph<DataType, DynamicEdges>;

}


#endif //RNAMAKE_NEW_NEW_GRAPH_H

















































