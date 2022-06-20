//
// Created by Joseph Yesselman on 1/15/18.
//

#include <data_structure/graph/graph_adjacency_list.h>

namespace data_structure::graph {

///////////////////////////////////////////////////////////////////////////////
// AdjList ////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// to increase readabiltiy
// D is the type of data stored in a node
// E is the type of the edge either FixedEdges or DynamicEdges

// Construction ///////////////////////////////////////////////////////////////
template <typename D, typename E>
AdjList<D, E>::AdjList(const AdjList &adj_list) {
  _copy_list(adj_list);
}

// Operators //////////////////////////////////////////////////////////////////
template <typename D, typename E>
AdjList<D, E> &AdjList<D, E>::operator=(const AdjList<D, E> &adj_list) {
  _copy_list(adj_list);
  return *this;
}

// Non const methods //////////////////////////////////////////////////////////
template <typename D, typename E>
void AdjList<D, E>::add_edge(const NodeIndexandEdge &nie1,
                             const NodeIndexandEdge &nie2) {
  if (std::is_same<E, DynamicEdges>::value) {
    _update_dyanmic_edges(nie1);
    _update_dyanmic_edges(nie2);
  }
  if (!edge_index_empty(nie1.ni, nie1.ei)) {
    String msg = "cannot add edge, " + nie1.get_str() +
                 " already is connected to -> " +
                 get_connected_node_info(nie1).to_str();
    base::log_and_throw<GraphException>(msg);
  }
  if (!edge_index_empty(nie2.ni, nie2.ei)) {
    throw GraphException("cannot add edge, " + nie2.get_str() +
                         " already is connected to -> " +
                         get_connected_node_info(nie2).to_str());
  }
  _edges[nie1.ni][nie1.ei] = new Connection{nie1.ni, nie2.ni, nie1.ei, nie2.ei};
  _edges[nie2.ni][nie2.ei] = new Connection{nie1.ni, nie2.ni, nie1.ei, nie2.ei};
}

// getters ////////////////////////////////////////////////////////////////////
template <typename D, typename E>
bool AdjList<D, E>::are_nodes_connected(Index n1, Index n2) const {
  if (_edges.find(n1) == _edges.end()) {
    return false;
  }
  auto const &edges = _edges.at(n1);
  for (auto const &edge : edges) {
    if (edge == nullptr) {
      continue;
    }
    if (edge->node_i == n2 or edge->node_j == n2) {
      return true;
    }
  }
  return false;
}

template <typename D, typename E>
NodeIndexandEdge
AdjList<D, E>::get_connected_node_info(const NodeIndexandEdge &nie) const {
  if (_edges.find(nie.ni) == _edges.end()) {
    String msg = "cannot find node of index: " + std::to_string(nie.ni);
    base::log_and_throw<GraphException>(msg);
  }
  if (_edges.find(nie.ni)->second.size() < nie.ei) {
    String msg =
        "node has fewer edges then requested one: " + std::to_string(nie.ni);
    base::log_and_throw<GraphException>(msg);
  }
  auto &edge = _edges.find(nie.ni)->second[nie.ei];
  auto partner_index = edge->get_partner(nie.ni);
  auto partner_end_index = edge->get_edge_index(partner_index);
  return NodeIndexandEdge{partner_index, partner_end_index};
}

// protected functions ////////////////////////////////////////////////////////
template <typename D, typename E>
void AdjList<D, E>::_copy_list(const AdjList<D, E> &adj_list) {
  _nodes = adj_list._nodes;
  _index = adj_list._index;
  for (auto const &kv : _nodes) {
    _edges[kv.first] = Connections(adj_list.get_node_edges(kv.first).size());
  }
  for (auto const &kv : adj_list._edges) {
    for (auto const &e : kv.second) {
      if (e == nullptr) {
        continue;
      }
      if (edge_between_nodes(e->node_i, e->node_j)) {
        continue;
      }
      add_edge(NodeIndexandEdge{e->node_i, e->edge_i},
               NodeIndexandEdge{e->node_j, e->edge_j});
    }
  }
}

} // namespace data_structure::graph
