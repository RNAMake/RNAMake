//
// Created by Joseph Yesselman on 1/15/18.
//

#ifndef RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H
#define RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H

#include <base/log.hpp>
#include <data_structure/graph/graph_base.h>

namespace data_structure::graph {

// NOTE this is a highly templated class, so it cannot be seperated into cpp
// adjecnecy list
// Data is the type of data stored in a node
// Edge is the type of the edge either FixedEdges or DynamicEdges
template <typename Data, typename Edge> class AdjacencyList {
public: // construction ///////////////////////////////////////////////////////
  inline AdjacencyList() = default;
  /// @brief copy constructor
  AdjacencyList(const AdjacencyList &adj_list) { _copy_list(adj_list); }
  /// @brief specifed destructor needs to make sure all pointers get deleted
  virtual ~AdjacencyList() {
    for (auto &kv : _connections) {
      for (auto &e : kv.second) {
        delete e;
      }
    }
  }

public: // operators //////////////////////////////////////////////////////////
  AdjacencyList &operator=(const AdjacencyList &adj_list) {
    _copy_list(adj_list);
    return *this;
  }

public: // iteration //////////////////////////////////////////////////////////
  typedef typename std::map<Index, Node<Data>>::const_iterator const_iterator;

  const_iterator begin() const { return _nodes.begin(); }
  const_iterator end() const { return _nodes.end(); }

public: // non const methods //////////////////////////////////////////////////
  void add_connection(const ConnectionPoint &cp1, const ConnectionPoint &cp2) {
    if (std::is_same<Edge, DynamicEdges>::value) {
      _update_dyanmic_edges(cp1);
      _update_dyanmic_edges(cp2);
    }
    if (!connection_point_empty(cp1.ni, cp1.ei)) {
      String msg = "cannot add edge, " + cp1.get_str() +
                   " already is connected to -> " +
                   get_connected_node_info(cp1).get_str();
      base::log_and_throw<GraphException>(msg);
    }
    if (!connection_point_empty(cp2.ni, cp2.ei)) {
      String msg = ("cannot add edge, " + cp2.get_str() +
                    " already is connected to -> " +
                    get_connected_node_info(cp2).get_str());
      base::log_and_throw<GraphException>(msg);
    }
    _connections[cp1.ni][cp1.ei] = new Connection{cp1, cp2};
    _connections[cp2.ni][cp2.ei] = new Connection{cp1, cp2};
  }

  virtual Index add_node(Data &d, Size n_edges) {
    _nodes.insert({_index, Node<Data>{d, _index}});
    _connections[_index] = Connections(n_edges);
    _index += 1;
    return _index - 1;
  }

  virtual void remove_connection(ConnectionPoint const &cp1,
                                 ConnectionPoint const &cp2) {
    if (connection_point_empty(cp1.ni, cp1.ei)) {
      String msg = "cannot remove edge, " + cp1.get_str() + " it is not filled";
      base::log_and_throw<GraphException>(msg);
    }
    if (connection_point_empty(cp2.ni, cp2.ei)) {
      String msg = "cannot remove edge, " + cp2.get_str() + " it is not filled";
      base::log_and_throw<GraphException>(msg);
    }
    // auto other_cp = get_connected_node_info(cp1);
    // if(other_cp != cp2) {

    //}
    delete _connections[cp1.ni][cp1.ei];
    delete _connections[cp2.ni][cp2.ei];
    _connections[cp1.ni][cp1.ei] = nullptr;
    _connections[cp2.ni][cp2.ei] = nullptr;
  }

  virtual void remove_node(Index ni) {
    auto edges = get_node_connections(ni);
    for (auto &e : edges) {
      if (e == nullptr) {
        continue;
      }
      Index partner_index = e->get_partner_index(ni);
      Index partner_end_index = e->get_edge_index(partner_index);
      delete _connections[partner_index][partner_end_index];
      _connections[partner_index][partner_end_index] = nullptr;
    }

    for (int i = 0; i < edges.size(); i++) {
      if (edges[i] == nullptr) {
        continue;
      }
      delete edges[i];
      _connections[ni][i] = nullptr;
    }
    _connections.erase(ni);
    _nodes.erase(ni);
    _update_next_index();
  }

public: // getters ////////////////////////////////////////////////////////////
  [[nodiscard]] bool are_nodes_connected(Index n1, Index n2) const {
    if (_connections.find(n1) == _connections.end()) {
      return false;
    }
    auto const &edges = _connections.at(n1);
    for (auto const &edge : edges) {
      if (edge == nullptr) {
        continue;
      }
      if(edge->cp1.ni == n1 and edge->cp2.ni == n2) {
        return true;
      }
      if(edge->cp1.ni == n2 and edge->cp2.ni == n1 ) {
        return true;
      }
    }
    return false;
  }
  /// @brief checks to see if an edge on a node is filled
  [[nodiscard]] bool connection_point_empty(Index ni, Index ei) const {
    if (_connections.find(ni) == _connections.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    if (_connections.find(ni)->second.size() < ei + 1) {
      String msg =
          "node has fewer edges then requested one: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    if (_connections.find(ni)->second[ei] == nullptr) {
      return true;
    } else {
      return false;
    }
  }
  /// @brief returns the node index and edge index of other node in connection
  [[nodiscard]] ConnectionPoint
  get_connected_node_info(const ConnectionPoint &nie) const {
    if (_connections.find(nie.ni) == _connections.end()) {
      String msg = "cannot find node of index: " + std::to_string(nie.ni);
      base::log_and_throw<GraphException>(msg);
    }
    if (_connections.find(nie.ni)->second.size() < nie.ei) {
      String msg =
          "node has fewer edges then requested one: " + std::to_string(nie.ni);
      base::log_and_throw<GraphException>(msg);
    }
    auto &edge = _connections.find(nie.ni)->second[nie.ei];
    auto partner_index = edge->get_partner_index(nie.ni);
    auto partner_end_index = edge->get_edge_index(partner_index);
    return ConnectionPoint{partner_index, partner_end_index};
  }
  /// @brief get all the connections a node has both filled and not
  [[nodiscard]] const Connections &get_node_connections(Index ni) const {
    if (_nodes.find(ni) == _nodes.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    return _connections.at(ni);
  }

  [[nodiscard]] Node<Data> const &get_node(Index ni) const {
    if (_nodes.find(ni) == _nodes.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    return _nodes.find(ni)->second;
  }

  [[nodiscard]] Data const &get_node_data(Index ni) const {
    if (_nodes.find(ni) == _nodes.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    return _nodes.find(ni)->second.get_data();
  }

  [[nodiscard]] size_t get_num_connections() const {
    auto n_edges = 0;
    for (auto const &kv : _connections) {
      for (auto const &e : kv.second) {
        if (e != nullptr) {
          n_edges += 1;
        }
      }
    }
    return (size_t)n_edges / 2;
  }

  [[nodiscard]] inline size_t get_num_nodes() const { return _nodes.size(); }

public: // setters ////////////////////////////////////////////////////////////
  void set_node_data(Index ni, Data &data) {
    if (_nodes.find(ni) == _nodes.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    _nodes.find(ni)->second.set_data(data);
  }

protected: // internal functions //////////////////////////////////////////////
  void _update_dyanmic_edges(ConnectionPoint const &nie) {
    if (_connections.find(nie.ni) == _connections.end()) {
      String msg = "cannot find node of index: " + std::to_string(nie.ni);
      base::log_and_throw<GraphException>(msg);
    }
    if (_connections.find(nie.ni)->second.size() < nie.ei + 1) {
      _connections.find(nie.ni)->second.resize(nie.ei + 1);
    }
  }

  void _update_next_index() {
    auto largest = 0;
    if (_nodes.size() == 0) {
      _index = 0;
      return;
    }
    for (auto const &kv : _nodes) {
      if (kv.first > largest) {
        largest = kv.first;
      }
    }
    _index = largest + 1;
  }

  void _copy_list(const AdjacencyList &adj_list) {
    _nodes = adj_list._nodes;
    _index = adj_list._index;
    for (auto const &kv : _nodes) {
      _connections[kv.first] =
          Connections(adj_list.get_node_connections(kv.first).size());
    }
    for (auto const &kv : adj_list._connections) {
      for (auto const &e : kv.second) {
        if (e == nullptr) {
          continue;
        }
        if (are_nodes_connected(e->cp1.ni, e->cp2.ni)) {
          continue;
        }
        add_connection(ConnectionPoint{e->cp1.ni, e->cp1.ei},
                       ConnectionPoint{e->cp2.ni, e->cp2.ei});
      }
    }
  }

protected: // members /////////////////////////////////////////////////////////
  /// @brief all the connections for each node, retrivable by node index
  std::map<Index, Connections> _connections = {};
  /// @brief all the nodes in the graph, retrivable by node index
  std::map<Index, Node<Data>> _nodes = {};
  /// @brief current index of next node
  Index _index = 0;
};

template <typename Data> using FixedEdged_AL = AdjacencyList<Data, FixedEdges>;

template <typename Data, typename Edge>
class DirectedAdjacencyList : public AdjacencyList<Data, Edge> {
public:
  typedef AdjacencyList<Data, Edge> BaseClass;

public:
  DirectedAdjacencyList() : BaseClass() {}

  DirectedAdjacencyList(const DirectedAdjacencyList &abj_list)
      : DirectedAdjacencyList() {
    this->parent_ = abj_list.parent_;
    this->_copy_list(abj_list);
  }

public:
  DirectedAdjacencyList &operator=(DirectedAdjacencyList const &abj_list) {
    this->parent_ = abj_list.parent_;
    this->_copy_list(abj_list);
    return *this;
  }

public:
  Index add_node(Data &d, Size n_edges) {
    auto ni = BaseClass::add_node(d, n_edges);
    parent_[ni] = -1;
    return ni;
  }

  Index add_node(Data &d, Size n_edges, Index n_end_index,
                 const ConnectionPoint &pie) {

    // not sure why this would happen but will catch anyway
    // expects<GraphException>(n_edges > n_end_index,
    //                        "n_edges must be greater than n_end_index");

    // parent should exist
    // expects<GraphException>(this->_nodes.find(pie.ni) != this->_nodes.end(),
    //                        "cannot add node to graph, parent with index: " +
    //                            std::to_string(pie.ni) + " does not exist");

    auto n_index = add_node(d, n_edges);
    BaseClass::add_connection(pie, ConnectionPoint{n_index, n_end_index});
    parent_[n_index] = pie.ni;
    return n_index;
  }

  void remove_node(Index ni) {
    BaseClass::remove_node(ni);
    for (auto &kv : parent_) {
      if (kv.second == ni) {
        kv.second = -1;
      }
    }
    parent_.erase(ni);
  }

public:
  [[nodiscard]] bool has_parent(Index n_index) const {
    if (parent_.at(n_index) == -1) {
      return false;
    } else {
      return true;
    }
  }

  [[nodiscard]] Index get_parent_index(Index n_index) const {

    // expects<GraphException>(
    //     parent_.at(n_index) != -1,
    //     "node does not have parent cannot get parent index");
    return parent_.find(n_index)->second;
  }

  [[nodiscard]] Index get_parent_end_index(Index n_index) const {
    auto parent_index = get_parent_index(n_index);
    auto &edges = this->get_node_edges(n_index);
    for (auto const &edge : edges) {
      if (edge == nullptr) {
        continue;
      }
      if (edge->partner(n_index) == parent_index) {
        return edge->end_index(parent_index);
      }
    }
    throw GraphException("cannot find parent end index!!!!");
  }

protected:
  std::map<Index, Index> parent_;
};

template <typename DataType>
using FixedEdged_DAL = DirectedAdjacencyList<DataType, FixedEdges>;

} // namespace data_structure::graph

#endif // RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H
