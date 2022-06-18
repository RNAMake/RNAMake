//
// Created by Joseph Yesselman on 1/15/18.
//

#ifndef RNAMAKE_NEW_GRAPH_BASE_H
#define RNAMAKE_NEW_GRAPH_BASE_H

#include <map>

#include <base/exception.hpp>
#include <base/log.hpp>
#include <base/types.hpp>

namespace data_structure::graph {

/*
 * Exception for graph
 */
class GraphException : public std::runtime_error {
public:
  /**
   * Standard constructor for ChainException
   * @param   message   Error message for chain
   */
  explicit GraphException(String const &message)
      : std::runtime_error(message) {}
};

// @brief holds data and the index of a given node
template <typename DataType> class Node {
public:
  inline Node(DataType &data, const Index index)
      : _data(std::move(data)), _index(index) {}

public:
  inline DataType const &get_data() const { return _data; }

  [[nodiscard]] inline Index get_index() const { return _index; }

public:
  inline void set_data(DataType & data) {
    _data = std::move(data);
  }

private:
  DataType _data;
  Index _index;
};

/// @brief each graph node has an index and N edges NodeIndexandEdge allows
/// specification of which node and edge position is involved in a connection
struct NodeIndexandEdge {
public:
  [[nodiscard]] String get_str() const {
    return "ni: " + std::to_string(ni) + " ei: " + std::to_string(ei);
  }

public:     // members
  Index ni; // node index
  Index ei; // edge index
};

struct NodeIndexandEdgeCompare {
  bool operator()(const NodeIndexandEdge &nie1,
                  const NodeIndexandEdge &nie2) const {
    if (nie1.ni != nie2.ni) {
      return nie1.ni > nie2.ni;
    } else {
      return nie1.ei > nie2.ei;
    }
  }
};

typedef std::map<NodeIndexandEdge, NodeIndexandEdge, NodeIndexandEdgeCompare>
    NodeIndexandEdgeMap;

struct Connection {
public:
  inline bool operator==(const Connection &other) const {
    return node_i == other.node_i && node_j == other.node_j &&
           edge_i == other.edge_i && edge_j == other.edge_j;
  }

public: // getters
  /// @brief get other node in the connection
  [[nodiscard]] Index get_partner(Index) const;
  /// @brief get index of edge a node is connected in 
  [[nodiscard]] Index get_edge_index(Index index) const;
  /// @brief
  [[nodiscard]] String get_str() const;

public:         // members
  Index node_i; // index of node 1 in graph
  Index node_j; // index of node 2 in graph
  Index edge_i; // index of edge of connection in node 1
  Index edge_j; // index of edge of connection in node 2
};

// typedefs
typedef std::vector<Connection const *> Connections;

// Edge Types
struct FixedEdges {};
struct DynamicEdges {};

} // namespace data_structure::graph

#endif // RNAMAKE_NEW_GRAPH_BASE_H
