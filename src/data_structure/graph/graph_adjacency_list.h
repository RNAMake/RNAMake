//
// Created by Joseph Yesselman on 1/15/18.
//

#ifndef RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H
#define RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H

#include <base/log.hpp>
#include <data_structure/graph/graph_base.h>

namespace data_structure::graph {

// adjecnecy list
// D is the type of data stored in a node
// E is the type of the edge either FixedEdges or DynamicEdges
template <typename D, typename E> class AdjList {
public: // construction ///////////////////////////////////////////////////////
  inline AdjList() = default;

  AdjList(const AdjList &);

  virtual ~AdjList() {
    for (auto &kv : _edges) {
      for (auto &e : kv.second) {
        // TODO can I just do "delete e"?
        if (e != nullptr) {
          delete e;
        }
      }
    }
  }

public: // operators //////////////////////////////////////////////////////////
  AdjList &operator=(const AdjList &);

public: // iteration //////////////////////////////////////////////////////////
  typedef typename std::map<Index, Node<D>>::const_iterator const_iterator;

  const_iterator begin() const { return _nodes.begin(); }
  const_iterator end() const { return _nodes.end(); }

public: // non const methods //////////////////////////////////////////////////
  void add_edge(const NodeIndexandEdge &, const NodeIndexandEdge &);

  virtual Index add_node(D &d, Size n_edges) {
    _nodes.insert({_index, Node<D>{d, _index}});
    _edges[_index] = Connections(n_edges);
    _index += 1;
    return _index - 1;
  }

  virtual void remove_edge(NodeIndexandEdge const &nie1,
                           NodeIndexandEdge const &nie2) {
    if (edge_index_empty(nie1.ni, nie1.ei)) {
      String msg =
          "cannot remove edge, " + nie1.get_str() + " it is not filled";
      base::log_and_throw<GraphException>(msg);
    }
    if (edge_index_empty(nie2.ni, nie2.ei)) {
      String msg =
          "cannot remove edge, " + nie2.get_str() + " it is not filled";
      base::log_and_throw<GraphException>(msg);
    }
    delete _edges[nie1.ni][nie1.ei];
    delete _edges[nie2.ni][nie2.ei];
    _edges[nie1.ni][nie1.ei] = nullptr;
    _edges[nie2.ni][nie2.ei] = nullptr;
  }

  virtual void remove_node(Index ni) {
    auto edges = get_node_connections(ni);
    for (auto &e : edges) {
      if (e == nullptr) {
        continue;
      }
      Index partner_index = e->get_partner(ni);
      Index partner_end_index = e->get_partner(partner_index);
      delete _edges[partner_index][partner_end_index];
      _edges[partner_index][partner_end_index] = nullptr;
    }

    for (int i = 0; i < edges.size(); i++) {
      if (edges[i] == nullptr) {
        continue;
      }
      delete edges[i];
      _edges[ni][i] = nullptr;
    }
    _edges.erase(ni);
    _nodes.erase(ni);
    _update_next_index();
  }

public: // getters ////////////////////////////////////////////////////////////
  [[nodiscard]] bool are_nodes_connected(Index, Index) const;
  /// @brief checks to see if an edge on a node is filled
  /// must be implemented in hpp as its called by virtual member
  [[nodiscard]] bool edge_index_empty(Index ni, Index ei) const {
    if (_edges.find(ni) == _edges.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    if (_edges.find(ni)->second.size() < ei) {
      String msg =
          "node has fewer edges then requested one: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    if (_edges.find(ni)->second[ei] == nullptr) {
      return true;
    } else {
      return false;
    }
  }

  [[nodiscard]] NodeIndexandEdge
  get_connected_node_info(const NodeIndexandEdge &) const;
  /// @brief get all the connections a node has both filled and not
  /// must be implemented in hpp as its called by virtual member
  [[nodiscard]] const Connections &get_node_connections(Index ni) const {
    if (_nodes.find(ni) == _nodes.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    return _edges.at(ni);
  }

  [[nodiscard]] Node<D> const &get_node(Index ni) const {
    if (_nodes.find(ni) == _nodes.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    return _nodes.find(ni)->second;
  }

  [[nodiscard]] D const &get_node_data(Index ni) const {
    if (_nodes.find(ni) == _nodes.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    return _nodes.find(ni)->second.get_data();
  }

  [[nodiscard]] size_t get_num_connections() const {
    auto n_edges = 0;
    for (auto const &kv : _edges) {
      for (auto const &e : kv.second) {
        if (e != nullptr) {
          n_edges += 1;
        }
      }
    }
    return (size_t)n_edges / 2;
  }

  [[nodiscard]] inline size_t get_num_nodes() const { return _nodes.size(); }

public: // setters
  void set_node_data(Index ni, D &data) {
    if (_nodes.find(ni) == _nodes.end()) {
      String msg = "cannot find node of index: " + std::to_string(ni);
      base::log_and_throw<GraphException>(msg);
    }
    _nodes.find(ni)->second.set_data(data);
  }

protected:
  void _update_dyanmic_edges(NodeIndexandEdge const &nie) {
    if (_edges.find(nie.ni) == _edges.end()) {
      String msg = "cannot find node of index: " + std::to_string(nie.ni);
      base::log_and_throw<GraphException>(msg);
    }
    if (_edges.find(nie.ni)->second.size() < nie.ei - 1) {
      _edges.find(nie.ni)->second.resize(nie.ei);
    }
  }

  void _update_next_index() {
    auto largest = 0;
    for (auto const &kv : _nodes) {
      if (kv.first > largest) {
        largest = kv.first;
      }
    }
    _index = largest + 1;
  }

  void _copy_list(const AdjList &);

protected:
  std::map<Index, Connections> _edges = {};
  std::map<Index, Node<D>> _nodes = {};
  Index _index = 0;
};

/*
template <typename DataType>
using FixedEdged_AL = AdjacencyList<DataType, FixedEdges>;

template <typename DataType, typename EdgeType>
class DirectedAdjacencyList : public AdjacencyList<DataType, EdgeType> {
public:
  typedef AdjacencyList<DataType, EdgeType> BaseClass;

public:
  DirectedAdjacencyList() : BaseClass() {}

  DirectedAdjacencyList(DirectedAdjacencyList const &abj_list)
      : DirectedAdjacencyList() {
    this->_nodes = abj_list._nodes;
    this->_index = abj_list._index;
    this->parent_ = abj_list.parent_;
    for (auto const &kv : this->_nodes) {
      this->_edges[kv.first] =
          std::vector<Edge const *>(abj_list.get_node_edges(kv.first).size());
    }

    for (auto const &kv : abj_list._edges) {
      for (auto const &e : kv.second) {
        if (e == nullptr) {
          continue;
        }
        if (this->edge_between_nodes(e->node_i, e->node_j)) {
          continue;
        }
        this->add_edge(NodeIndexandEdge{e->node_i, e->edge_i},
                       NodeIndexandEdge{e->node_j, e->edge_j});
      }
    }
  }

public:
  DirectedAdjacencyList &operator=(DirectedAdjacencyList const &abj_list) {
    this->_nodes = abj_list._nodes;
    this->_index = abj_list._index;
    this->parent_ = abj_list.parent_;
    for (auto const &kv : this->_nodes) {
      this->_edges[kv.first] =
          std::vector<Edge const *>(abj_list.get_node_edges(kv.first).size());
    }

    for (auto const &kv : abj_list._edges) {
      for (auto const &e : kv.second) {
        if (e == nullptr) {
          continue;
        }
        if (this->edge_between_nodes(e->node_i, e->node_j)) {
          continue;
        }
        this->add_edge(NodeIndexandEdge{e->node_i, e->edge_i},
                       NodeIndexandEdge{e->node_j, e->edge_j});
      }
    }
    return *this;
  }

public:
  Index add_node(DataType const &d, Size n_edges) {
    auto ni = BaseClass::add_node(d, n_edges);
    parent_[ni] = -1;
    return ni;
  }

  Index add_node(DataType const &d, Size n_edges, Index n_end_index,
                 NodeIndexandEdge const &pie) {

    // not sure why this would happen but will catch anyway
    expects<GraphException>(n_edges > n_end_index,
                            "n_edges must be greater than n_end_index");

    // parent should exist
    expects<GraphException>(
        this->_nodes.find(pie.ni) != this->_nodes.end(),
        "cannot add node to graph, parent with index: " +
            std::to_string(pie.ni) + " does not exist");

    auto n_index = add_node(d, n_edges);
    BaseClass::add_edge(pie, NodeIndexandEdge{n_index, n_end_index});
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
  bool has_parent(Index n_index) const {
    if (parent_.at(n_index) == -1) {
      return false;
    } else {
      return true;
    }
  }

  Index get_parent_index(Index n_index) const {

    expects<GraphException>(
        parent_.at(n_index) != -1,
        "node does not have parent cannot get parent index");
    return parent_.find(n_index)->second;
  }

  Index get_parent_end_index(Index n_index) const {
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
 */

} // namespace data_structure::graph

#endif // RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H
