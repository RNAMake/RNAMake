//
// Created by Joseph Yesselman on 1/15/18.
//

#include <data_structure/graph/graph_base.h>

namespace data_structure::graph {

// Connection functions ///////////////////////////////////////////////////////
Index Connection::get_partner(Index index) const {
  if (index == node_i) {
    return node_j;
  } else if (index == node_j) {
    return node_i;
  } else {
    String msg = "cannot call partner if not a node in this edge";
    base::log_and_throw<GraphException>(msg);
    return -1; // stop warning
  }
}

Index Connection::get_edge_index(Index index) const {
  if (index == node_i) {
    return edge_i;
  } else if (index == node_j) {
    return edge_j;
  } else {
    String msg = "cannot call end_index if not a node in this edge";
    base::log_and_throw<GraphException>(msg);
    return -1; // stop warning
  }
}

String Connection::get_str() const {
  return "ni: " + std::to_string(node_i) + " nj: " + std::to_string(node_j) +
         " ei: " + std::to_string(edge_i) + " ej: " + std::to_string(edge_j);
}

} // namespace data_structure::graph