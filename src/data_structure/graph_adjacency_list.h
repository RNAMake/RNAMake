//
// Created by Joseph Yesselman on 1/15/18.
//

#ifndef RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H
#define RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H

#include <base/log.hpp>
#include <data_structure/graph_base.h>

namespace data_structure {

template <typename DataType, typename EdgeType> class AdjacencyList {
public:
  AdjacencyList()
      : edges_(std::map<Index, std::vector<Edge const *>>()),
        nodes_(std::map<Index, Node<DataType>>()), index_(0) {}

  AdjacencyList(AdjacencyList const &abj_list) {
    nodes_ = abj_list.nodes_;
    index_ = abj_list.index_;
    for (auto const &kv : nodes_) {
      edges_[kv.first] =
          std::vector<Edge const *>(abj_list.get_node_edges(kv.first).size());
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

  virtual ~AdjacencyList() {
    for (auto &kv : edges_) {
      for (auto &e : kv.second) {
        if (e != nullptr) {
          delete e;
        }
      }
    }
  }

public:
  AdjacencyList &operator=(AdjacencyList const &abj_list) {
    this->nodes_ = abj_list.nodes_;
    this->index_ = abj_list.index_;
    for (auto const &kv : nodes_) {
      this->edges_[kv.first] =
          std::vector<Edge const *>(abj_list.get_node_edges(kv.first).size());
    }
    for (auto const &kv : abj_list.edges_) {
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
  typedef
      typename std::map<Index, Node<DataType>>::const_iterator const_iterator;

  const_iterator begin() const { return nodes_.begin(); }
  const_iterator end() const { return nodes_.end(); }

public:
  virtual Index add_node(DataType const &d, Size n_edges) {
    nodes_.insert(
        std::pair<int, Node<DataType>>(index_, Node<DataType>(d, index_)));
    edges_[index_] = std::vector<Edge const *>(n_edges);
    index_ += 1;
    return index_ - 1;
  }

  virtual void add_edge(NodeIndexandEdge const &nie1,
                        NodeIndexandEdge const &nie2) {

    if (std::is_same<EdgeType, DynamicEdges>::value) {
      _update_dyanmic_edges(nie1);
      _update_dyanmic_edges(nie2);
    }

    if (!edge_index_empty(nie1.node_index, nie1.edge_index)) {
      throw GraphException("cannot add edge, " + nie1.to_str() +
                           " already is connected to -> " +
                           get_connected_node_info(nie1).to_str());
    }
    if (!edge_index_empty(nie2.node_index, nie2.edge_index)) {
      throw GraphException("cannot add edge, " + nie2.to_str() +
                           " already is connected to -> " +
                           get_connected_node_info(nie2).to_str());
    }

    edges_[nie1.node_index][nie1.edge_index] = new Edge(
        nie1.node_index, nie2.node_index, nie1.edge_index, nie2.edge_index);
    edges_[nie2.node_index][nie2.edge_index] = new Edge(
        nie1.node_index, nie2.node_index, nie1.edge_index, nie2.edge_index);
  }

  virtual void remove_node(Index ni) {
    auto edges = get_node_edges(ni);
    for (auto &e : edges) {
      if (e == nullptr) {
        continue;
      }
      auto partner_index = e->partner(ni);
      auto partner_end_index = e->end_index(partner_index);
      delete edges_[partner_index][partner_end_index];
      edges_[partner_index][partner_end_index] = nullptr;
    }

    for (int i = 0; i < edges.size(); i++) {
      if (edges[i] == nullptr) {
        continue;
      }
      delete edges[i];
      edges_[ni][i] = nullptr;
    }
    edges_.erase(ni);
    nodes_.erase(ni);
    _update_next_index();
  }

  virtual void remove_edge(NodeIndexandEdge const &nie1,
                           NodeIndexandEdge const &nie2) {

    if (edge_index_empty(nie1.node_index, nie1.edge_index)) {
      throw GraphException("cannot remove edge, " + nie1.to_str() +
                           " it is not filled");
    }
    if (edge_index_empty(nie2.node_index, nie2.edge_index)) {
      throw GraphException("cannot remove edge, " + nie2.to_str() +
                           " it is not filled");
    }

    delete edges_[nie1.node_index][nie1.edge_index];
    delete edges_[nie2.node_index][nie2.edge_index];

    edges_[nie1.node_index][nie1.edge_index] = nullptr;
    edges_[nie2.node_index][nie2.edge_index] = nullptr;
  }

public:
  inline size_t get_num_nodes() const { return nodes_.size(); }

  size_t get_num_edges() const {
    auto n_edges = 0;
    for (auto const &kv : edges_) {
      for (auto const &e : kv.second) {
        if (e != nullptr) {
          n_edges += 1;
        }
      }
    }
    return (size_t)n_edges / 2;
  }

  std::vector<Edge const *> const &get_node_edges(Index ni) const {
    expects<GraphException>(edges_.find(ni) != edges_.end(),
                            "this node has no edges : " + std::to_string(ni));
    return edges_.at(ni);
  }

  DataType const &get_node_data(Index ni) const {
    std::cout << "I am in get node data" << std::endl;
    expects<GraphException>(nodes_.find(ni) != nodes_.end(),
                            "cannot find node of index: " + std::to_string(ni));
    return nodes_.find(ni)->second.data();
  }

  DataType &get_node_data(Index ni) {
    std::cout << "I am in get node data" << std::endl;
    expects<GraphException>(nodes_.find(ni) != nodes_.end(),
                            "cannot find node of index: " + std::to_string(ni));
    return nodes_.find(ni)->second.data();
  }

  Node<DataType> const &get_node(Index ni) const {
    expects<GraphException>(nodes_.find(ni) != nodes_.end(),
                            "cannot find node of index: " + std::to_string(ni));
    return nodes_.find(ni)->second;
  }

  Node<DataType> &get_node(Index ni) {
    expects<GraphException>(nodes_.find(ni) != nodes_.end(),
                            "cannot find node of index: " + std::to_string(ni));
    return nodes_.find(ni)->second;
  }

  NodeIndexandEdge get_connected_node_info(NodeIndexandEdge const &nei) const {

    expects<GraphException>(edges_.find(nei.node_index) != edges_.end() ||
                                edges_.find(nei.node_index)->second.size() >
                                    nei.edge_index - 1,
                            "node has fewer edges then requested one");

    auto edge = edges_.find(nei.node_index)->second[nei.edge_index];
    auto partner_index = edge->partner(nei.node_index);
    auto partner_end_index = edge->end_index(partner_index);
    return NodeIndexandEdge{partner_index, partner_end_index};
  }

  bool edge_between_nodes(Index n1, Index n2) const {
    if (edges_.find(n1) == edges_.end()) {
      return false;
    }
    auto &edges = edges_.at(n1);
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

public:
  bool edge_index_empty(Index ni, Index ei) const {

    expects<GraphException>(edges_.find(ni) != edges_.end() ||
                                edges_.find(ni)->second.size() > ei - 1,
                            "node has fewer edges then requested one");

    if (edges_.find(ni)->second[ei] == nullptr) {
      return true;
    } else {
      return false;
    }
  }

protected:
  void _update_dyanmic_edges(NodeIndexandEdge const &nie) {
    expects<GraphException>(
        edges_.find(nie.node_index) != edges_.end(),
        "node with index: " + std::to_string(nie.node_index) +
            " does not exist!");

    if (edges_.find(nie.node_index)->second.size() < nie.edge_index - 1) {
      edges_.find(nie.node_index)->second.resize(nie.edge_index);
    }
  }

  void _update_next_index() {
    auto largest = 0;
    for (auto const &kv : nodes_) {
      if (kv.first > largest) {
        largest = kv.first;
      }
    }
    index_ = largest + 1;
  }

protected:
  std::map<Index, std::vector<Edge const *>> edges_;
  std::map<Index, Node<DataType>> nodes_;
  Index index_;
};

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
    this->nodes_ = abj_list.nodes_;
    this->index_ = abj_list.index_;
    this->parent_ = abj_list.parent_;
    for (auto const &kv : this->nodes_) {
      this->edges_[kv.first] =
          std::vector<Edge const *>(abj_list.get_node_edges(kv.first).size());
    }

    for (auto const &kv : abj_list.edges_) {
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
    this->nodes_ = abj_list.nodes_;
    this->index_ = abj_list.index_;
    this->parent_ = abj_list.parent_;
    for (auto const &kv : this->nodes_) {
      this->edges_[kv.first] =
          std::vector<Edge const *>(abj_list.get_node_edges(kv.first).size());
    }

    for (auto const &kv : abj_list.edges_) {
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
        this->nodes_.find(pie.node_index) != this->nodes_.end(),
        "cannot add node to graph, parent with index: " +
            std::to_string(pie.node_index) + " does not exist");

    auto n_index = add_node(d, n_edges);
    BaseClass::add_edge(pie, NodeIndexandEdge{n_index, n_end_index});
    parent_[n_index] = pie.node_index;
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

} // namespace data_structure

#endif // RNAMAKE_NEW_GRAPH_ADJACENCY_LIST_H
