//
// Created by Joseph Yesselman on 1/15/18.
//

#include <data_structure/graph/graph_adjacency_list.h>

namespace data_structure::graph {

// Construction ///////////////////////////////////////////////////////////////
template <typename DataType, typename EdgeType>
AdjacencyList<DataType, EdgeType>::AdjacencyList(
    const AdjacencyList &abj_list) {
  nodes_ = abj_list.nodes_;
  index_ = abj_list.index_;
  for (auto const &kv : nodes_) {
    edges_[kv.first] = Connections(abj_list.get_node_edges(kv.first).size());
  }
  for (auto const &kv : abj_list.edges_) {
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
