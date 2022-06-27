//
// Created by Joseph Yesselman on 1/15/18.
//

#include <data_structure/graph/graph_base.h>

namespace data_structure::graph {

// Connection functions ///////////////////////////////////////////////////////
const ConnectionPoint &
Connection::get_partner(const ConnectionPoint &cp) const {
  if (cp == cp1) {
    return cp2;
  } else if (cp == cp2) {
    return cp1;
  } else {
    String msg = "cannot call partner if not a node in this connection";
    base::log_and_throw<GraphException>(msg);
    return cp1; // stop warning
  }
}

Index Connection::get_partner_index(Index index) const {
  if (index == cp1.ni) {
    return cp2.ni;
  } else if (index == cp2.ni) {
    return cp1.ni;
  } else {
    String msg = "cannot call partner if not a node in this connection";
    base::log_and_throw<GraphException>(msg);
    return -1; // stop warning
  }
}

Index Connection::get_edge_index(Index index) const {
  if (index == cp1.ni) {
    return cp1.ei;
  } else if (index == cp2.ni) {
    return cp2.ei;
  } else {
    String msg = "cannot call end_index if not a node in this connection";
    base::log_and_throw<GraphException>(msg);
    return -1; // stop warning
  }
}

String Connection::get_str() const {
  return "ni: " + std::to_string(cp1.ni) + " nj: " + std::to_string(cp1.ei) +
         " ei: " + std::to_string(cp2.ni) + " ej: " + std::to_string(cp2.ei);
}

} // namespace data_structure::graph