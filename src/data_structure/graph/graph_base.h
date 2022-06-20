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
  inline Node(DataType &data, Index index)
      : _data(std::move(data)), _index(index) {}

public:
  inline DataType const &get_data() const { return _data; }

  [[nodiscard]] inline Index get_index() const { return _index; }

public:
  inline void set_data(DataType &data) { _data = std::move(data); }

private:
  DataType _data;
  Index _index;
};

/// @brief each graph node has an index and N edges ConnectionPoint allows
/// specification of which node and edge position is involved in a connection
struct ConnectionPoint {
public:
  inline bool operator==(const ConnectionPoint &other) const {
    return ni == other.ni && ei == other.ei;
  }

  inline bool operator!=(const ConnectionPoint &other) const {
    return !(ni == other.ni && ei == other.ei);
  }

public:
  [[nodiscard]] String get_str() const {
    return "ni: " + std::to_string(ni) + " ei: " + std::to_string(ei);
  }

public:     // members
  Index ni; // node index
  Index ei; // edge index
};

struct ConnectionPointCompare {
  bool operator()(const ConnectionPoint &nie1,
                  const ConnectionPoint &nie2) const {
    if (nie1.ni != nie2.ni) {
      return nie1.ni > nie2.ni;
    } else {
      return nie1.ei > nie2.ei;
    }
  }
};

typedef std::map<ConnectionPoint, ConnectionPoint, ConnectionPointCompare>
    ConnectionPointMap;

struct Connection {
public:
  inline bool operator==(const Connection &other) const {
    return cp1 == other.cp1 && cp2 == other.cp2;
  }

public: // getters
  /// @brief get other connection point
  [[nodiscard]] const ConnectionPoint &
  get_partner(const ConnectionPoint &) const;
  /// @brief get other node in the connection
  [[nodiscard]] Index get_partner_index(Index) const;
  /// @brief get index of edge a node is connected in
  [[nodiscard]] Index get_edge_index(Index index) const;
  /// @brief
  [[nodiscard]] String get_str() const;

public: // members
  ConnectionPoint cp1;
  ConnectionPoint cp2;
};

// typedefs
typedef std::vector<Connection const *> Connections;

// Edge Types
struct FixedEdges {};
struct DynamicEdges {};

} // namespace data_structure::graph

#endif // RNAMAKE_NEW_GRAPH_BASE_H
