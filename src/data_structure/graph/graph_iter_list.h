//
// Created by Joseph Yesselman on 1/19/18.
//

#ifndef RNAMAKE_NEW_GRAPH_ITER_LIST_H
#define RNAMAKE_NEW_GRAPH_ITER_LIST_H

// standard headers
#include <queue>

// RNAMake headers
#include <base/types.hpp>
#include <data_structure/graph/graph_adjacency_list.h>
#include <data_structure/graph/graph_base.h>

namespace data_structure::graph {

template <typename Data, typename AdjacencyList> class IterList {
public:
  struct VisitedNode {
  public:
    inline VisitedNode(std::shared_ptr<VisitedNode> nparent, Index nindex)
        : parent(nparent), index(nindex) {}

  public:
    bool index_in_path(Index i) {
      if (index == i) {
        return true;
      }
      auto current = parent;
      while (current != nullptr) {
        if (current->index == i) {
          return true;
        }
        current = current->parent;
      }
      return false;
    }

    int path_length() {
      auto path_length = 1;
      auto current = parent;
      while (current != nullptr) {
        path_length += 1;
        current = current->parent;
      }
      return path_length;
    }

  public:
    std::shared_ptr<VisitedNode> parent;
    Index index;
  };

  typedef std::shared_ptr<VisitedNode> VisitedNodeOP;

public:
  IterList() = default;

  virtual ~IterList() = default;

public:
  typedef
      typename std::vector<Node<Data> const *>::const_iterator const_iterator;

  const_iterator begin() const noexcept { return _iter_list.begin(); }
  const_iterator end() const noexcept { return _iter_list.end(); }

public:
  virtual void transversal(AdjacencyList &adj_list, Index start_n) {
    _iter_list.resize(0);
    _seen = std::map<Index, int>();
    auto neighbors = std::vector<Index>();
    neighbors.reserve(10);
    _open.push(start_n);
    _seen[start_n] = 1;
    while (!_open.empty()) {
      auto current = _open.front();
      _open.pop();
      const Node<Data> *p = &adj_list.get_node(current);
      _iter_list.push_back(p);
      _get_neighbors(current, adj_list, neighbors);
      for (auto const &n : neighbors) {
        _open.push(n);
      }
    }
    // ensure you find all nodes not connected or in another sub graph
    for (auto const &kv : adj_list) {
      if (_seen.find(kv.first) != _seen.end()) {
        continue;
      }
      auto new_start = kv.first;
      _open.push(new_start);
      _seen[new_start] = 1;
      while (!_open.empty()) {
        auto current = _open.front();
        _open.pop();
        const Node<Data> *p = &adj_list.get_node(current);
        _iter_list.push_back(p);
        _get_neighbors(current, adj_list, neighbors);
        for (auto const &n : neighbors) {
          _open.push(n);
        }
      }
    }
  }

  virtual void path_transversal(AdjacencyList &adj_list, Index start_n,
                                Index end_n) {
    _iter_list.resize(0);
    _seen = std::map<Index, int>();
    auto neighbors = std::vector<Index>();
    neighbors.reserve(10);
    auto start_vn = std::make_shared<VisitedNode>(nullptr, start_n);
    auto current_round = std::vector<VisitedNodeOP>{start_vn};
    auto next_round = std::vector<VisitedNodeOP>();
    auto end_vn = VisitedNodeOP(nullptr);
    while (current_round.size() > 0) {
      for (auto vn : current_round) {
        _get_neighbors_path(vn, adj_list, neighbors);
        for (auto const &ni : neighbors) {
          auto new_vn = std::make_shared<VisitedNode>(vn, ni);
          if (ni == end_n) {
            if (end_vn == nullptr) {
              end_vn = new_vn;
            } else if (end_vn->path_length() > new_vn->path_length()) {
              end_vn = new_vn;
            }
          } else {
            next_round.push_back(new_vn);
          }
        }
      }
      current_round = next_round;
      next_round = std::vector<VisitedNodeOP>();
    }
    if (end_vn == nullptr) {
      throw GraphException(
          "there is no path between nodes: " + std::to_string(start_n) +
          " and " + std::to_string(end_n) + " !");
    }
    auto path_back = std::vector<Index>();
    while (end_vn != nullptr) {
      path_back.push_back(end_vn->index);
      end_vn = end_vn->parent;
    }
    int pos = path_back.size() - 1;
    for (int i = pos; i >= 0; i--) {
      _iter_list.push_back(&adj_list.get_node(path_back[i]));
    }
  }

protected:
  virtual void _get_neighbors(Index ni, AdjacencyList &adj_list,
                              std::vector<Index> &neighbors) {
    neighbors.resize(0);
    auto &connections = adj_list.get_node_connections(ni);
    for (auto const &c : connections) {
      if (c == nullptr) {
        continue;
      }
      if (_seen.find(c->get_partner_index(ni)) != _seen.end()) {
        continue;
      }
      neighbors.push_back(c->get_partner_index(ni));
      _seen[c->get_partner_index(ni)];
    }
  }

  virtual void _get_neighbors_path(VisitedNodeOP vn, AdjacencyList &adj_list,
                                   std::vector<Index> &neighbors) {
    neighbors.resize(0);
    auto &connections = adj_list.get_node_connections(vn->index);
    for (auto const &c : connections) {
      if (c == nullptr) {
        continue;
      }
      if (vn->index_in_path(c->get_partner_index(vn->index))) {
        continue;
      }
      neighbors.push_back(c->get_partner_index(vn->index));
    }
  }

protected:
  std::vector<const Node<Data> *> _iter_list = {};
  std::queue<Index> _open = {};
  std::map<Index, int> _seen = {};
};

template <typename Data, typename AdjacencyList>
class DirectedIterList : public IterList<Data, AdjacencyList> {
public:
  typedef IterList<Data, AdjacencyList> BaseClass;
  typedef typename BaseClass::VisitedNodeOP VisitedNodeOP;

public:
  DirectedIterList() : BaseClass() {}

public:
  virtual void transversal(AdjacencyList &adj_list, Index start_n) {
    this->_iter_list.resize(0);

    this->_seen = std::map<Index, int>();
    auto neighbors = std::vector<Index>();
    neighbors.reserve(10);

    this->_open.push(start_n);
    this->_seen[start_n] = 1;

    while (this->_open.size() > 0) {
      auto current = this->_open.front();
      this->_open.pop();
      const Node<Data> *p = &adj_list.get_node(current);
      this->_iter_list.push_back(p);
      _get_neighbors(current, adj_list, neighbors);
      for (auto const &n : neighbors) {
        this->_open.push(n);
      }
    }

    // ensure you find all nodes not connected or in another sub graph
    bool found = true;
    while (found) {
      found = false;
      for (auto const &kv : adj_list) {
        if (this->_seen.find(kv.first) != this->_seen.end()) {
          continue;
        }
        if (adj_list.has_parent(kv.first)) {
          continue;
        }
        found = true;
        auto new_start = kv.first;
        this->_open.push(new_start);
        this->_seen[new_start] = 1;

        while (this->_open.size() > 0) {
          auto current = this->_open.front();
          this->_open.pop();
          const Node<Data> *p = &adj_list.get_node(current);
          this->_iter_list.push_back(p);
          _get_neighbors(current, adj_list, neighbors);
          for (auto const &n : neighbors) {
            this->_open.push(n);
          }
        }
      }
    }
  }

  void sub_graph_transversal(AdjacencyList &adj_list, Index start_n,
                             Index end_n) {
    BaseClass::path_transversal(adj_list, start_n, end_n);
    auto neighbors = std::vector<Index>();
    neighbors.reserve(10);
    int pos = 0;
    while (pos < this->_iter_list.size()) {
      _get_neighbors(this->_iter_list[pos]->get_index(), adj_list, neighbors);
      for (auto const &ni : neighbors) {
        auto found = false;
        for (auto const &n : this->_iter_list) {
          if (n->get_index() == ni) {
            found = true;
            break;
          }
        }
        if (!found) {
          this->_iter_list.push_back(&adj_list.get_node(ni));
        }
      }

      pos += 1;
    }
  }

protected:
  virtual void _get_neighbors(Index ni, AdjacencyList &adj_list,
                              std::vector<Index> &neighbors) {
    neighbors.resize(0);
    auto &connections = adj_list.get_node_connections(ni);
    for (auto const &c : connections) {
      if (c == nullptr) {
        continue;
      }
      auto pi = c->get_partner_index(ni);
      if (adj_list.has_parent(pi) && adj_list.get_parent_index(pi) == ni) {
        if (this->_seen.find(c->get_partner_index(ni)) != this->_seen.end()) {
          continue;
        }
        neighbors.push_back(c->get_partner_index(ni));
        this->_seen[c->get_partner_index(ni)];
      }
    }
  }

  virtual void _get_neighbors_path(VisitedNodeOP vn, AdjacencyList &adj_list,
                                   std::vector<Index> &neighbors) {
    neighbors.resize(0);
    auto &connections = adj_list.get_node_connections(vn->index);
    for (auto const &c : connections) {
      if (c == nullptr) {
        continue;
      }
      auto pi = c->get_partner_index(vn->index);
      if (adj_list.has_parent(pi) &&
          adj_list.get_parent_index(pi) == vn->index) {
        if (vn->index_in_path(c->get_partner_index(vn->index))) {
          continue;
        }
        neighbors.push_back(c->get_partner_index(vn->index));
      }
    }
  }
};

} // namespace data_structure::graph

#endif // RNAMAKE_NEW_GRAPH_ITER_LIST_H
