//
// Created by Joseph Yesselman on 1/14/18.
//

#ifndef RNAMAKE_NEW_NEW_GRAPH_H
#define RNAMAKE_NEW_NEW_GRAPH_H

#include <queue>

#include <base/types.hpp>
#include <data_structure/graph/graph_adjacency_list.h>
#include <data_structure/graph/graph_base.h>
#include <data_structure/graph/graph_iter_list.h>

namespace data_structure::graph {

template <typename Data, typename AdjacencyList, typename IterList>
class _Graph {
public: // construction ///////////////////////////////////////////////////////
  _Graph() = default;

  virtual ~_Graph() = default;

public: // iteration //////////////////////////////////////////////////////////
  typedef typename IterList::const_iterator const_iterator;

  const_iterator begin() const noexcept { return iter_list_.begin(); }
  const_iterator end() const noexcept { return iter_list_.end(); }

public: // operators //////////////////////////////////////////////////////////
  Data &operator[](Index ni) { return get_node_data(ni); }

public: // transversal ////////////////////////////////////////////////////////
  void setup_transversal(Index start_n) {
    iter_list_.transversal(adjacency_list_, start_n);
  }

  void setup_path_transversal(Index start_n, Index end_n) {
    iter_list_.path_transversal(adjacency_list_, start_n, end_n);
  }

public: // node and connection management  ////////////////////////////////////
  virtual inline Index add_node(Data &d, Size n_edges) {
    return adjacency_list_.add_node(d, n_edges);
  }

  inline void add_connection(const ConnectionPoint &nie1,
                             const ConnectionPoint &nie2) {
    return adjacency_list_.add_connection(nie1, nie2);
  }

  inline void remove_node(Index ni) { return adjacency_list_.remove_node(ni); }

  inline void remove_connection(const ConnectionPoint &nie1,
                                const ConnectionPoint &nie2) {
    return adjacency_list_.remove_connection(nie1, nie2);
  }

public:
  inline size_t get_num_nodes() const {
    return adjacency_list_.get_num_nodes();
  }

  inline size_t get_num_connections() const {
    return adjacency_list_.get_num_connections();
  }

  inline const Connections &get_node_connections(Index ni) const {
    return adjacency_list_.get_node_connections(ni);
  }

  inline const Node<Data> &get_node(Index ni) const {
    return adjacency_list_.get_node(ni);
  }

  inline const Data &get_node_data(Index ni) const {
    return adjacency_list_.get_node_data(ni);
  }

  inline Data &get_node_data(Index ni) {
    return adjacency_list_.get_node_data(ni);
  }

  inline ConnectionPoint
  get_paired_connection_point(ConnectionPoint const &nei) const {
    return adjacency_list_.get_paired_connection_point(nei);
  }

  inline bool are_nodes_connected(Index n1, Index n2) const {
    return adjacency_list_.are_nodes_connected(n1, n2);
  }

  inline bool connection_point_empty(Index ni, Index ei) const {
    return adjacency_list_.connection_point_empty(ni, ei);
  }

public: // setters ////////////////////////////////////////////////////////////
  inline void set_node_data(Index ni, Data &d) {
    adjacency_list_.set_node_data(ni, d);
  }

protected:
  AdjacencyList adjacency_list_;
  mutable IterList iter_list_; // needs to update to iterate
};

template <typename Data, typename Edge>
using _UndirectedGraph = _Graph<Data, AdjacencyList<Data, Edge>,
                                IterList<Data, AdjacencyList<Data, Edge>>>;

template <typename Data, typename Edge>
using _DirectedGraph =
    _Graph<Data, DirectedAdjacencyList<Data, Edge>,
           DirectedIterList<Data, DirectedAdjacencyList<Data, Edge>>>;

template <typename DataType, typename EdgeType>
class UndirectedGraph : public _UndirectedGraph<DataType, EdgeType> {
public:
  UndirectedGraph() = default;

  UndirectedGraph(UndirectedGraph const &g) {
    this->adjacency_list_ = g.adjacency_list_;
  }
};

template <typename DataType>
using FixedEdgeUndirectedGraph = UndirectedGraph<DataType, FixedEdges>;

template <typename DataType>
using DynamicEdgedUndirectedGraph = UndirectedGraph<DataType, DynamicEdges>;

template <typename DataType, typename EdgeType>
class DirectedGraph : public _DirectedGraph<DataType, EdgeType> {
public:
  typedef _DirectedGraph<DataType, EdgeType> BaseClass;

public:
  DirectedGraph() = default;

  DirectedGraph(const DirectedGraph &g) : BaseClass() {
    this->adjacency_list_ = g.adjacency_list_;
  }

public:
  inline void setup_sub_graph_transversal(Index start_n, Index end_n) {
    return this->iter_list_.sub_graph_transversal(this->adjacency_list_,
                                                  start_n, end_n);
  }

public:
  inline Index add_node(DataType &d, Size n_edges, Index n_end_index,
                        const ConnectionPoint &pie) {
    return this->adjacency_list_.add_node(d, n_edges, n_end_index, pie);
  }

  inline Index add_node(DataType &d, Size n_edges) {
    return this->adjacency_list_.add_node(d, n_edges);
  }

public:
  [[nodiscard]] inline bool has_parent(Index ni) const {
    return this->adjacency_list_.has_parent(ni);
  }

  [[nodiscard]] inline Index get_parent_index(Index ni) const {
    return this->adjacency_list_.get_parent_index(ni);
  }

  [[nodiscard]] inline Index get_parent_end_index(Index ni) const {
    return this->adjacency_list_.get_parent_end_index(ni);
  }

public:
  Indexes get_root_indexes() {
    auto roots = Indexes();
    for (auto const &kv : this->adjacency_list_) {
      if (!this->adjacency_list_.has_parent(kv.first)) {
        roots.push_back(kv.first);
      }
    }
    return roots;
  }
};

template <typename DataType>
using FixedEdgeDirectedGraph = DirectedGraph<DataType, FixedEdges>;

template <typename DataType>
using DynamicEdgedDirectedGraph = DirectedGraph<DataType, DynamicEdges>;

} // namespace data_structure::graph

#endif // RNAMAKE_NEW_NEW_GRAPH_H
