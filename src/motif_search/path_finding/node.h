//
// Created by Joseph Yesselman on 3/18/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_NODE_H
#define RNAMAKE_NEW_PATH_FINDING_NODE_H

#include <queue>

#include <motif/motif_state.h>

namespace motif_search {
namespace path_finding {

class Node;
typedef std::shared_ptr<Node> NodeOP;

class Node {
public:
  inline Node(motif::MotifStateOP ms, NodeOP parent, float score, int level,
              int parent_end_index, int node_type)
      : _state(ms), _parent(parent), _score(score), _level(level),
        _parent_end_index(parent_end_index), _node_type(node_type),
        _size(ms->size()), _ss_score(ms->score()) {
    if (_parent != nullptr) {
      _ss_score += _parent->ss_score();
      // -2 for shared base pair
      _size += parent->size() - 2;
    }
  }

  ~Node() {}

public: // getters
  int get_level() const { return _level; }

  inline int get_size() const { return _size; }

  inline int get_node_type() const { return _node_type; }

  inline float get_ss_score() const { return _ss_score; }

  inline float get_score() const { return _score; }

  inline int get_parent_end_index() const { return _parent_end_index; }

  inline motif::MotifStateOP get_state() const { return _state; }

  inline NodeOP get_parent() const { return _parent; }

private:
  NodeOP _parent;
  motif::MotifStateOP _state;
  int _parent_end_index, _level, _size, _node_type;
  float _ss_score, _score;
};

struct NodeCompare {
  bool operator()(NodeOP node1, NodeOP node2) {

    if (node1->score() > node2->score()) {
      return true;
    } else {
      return false;
    }
  }
};

typedef std::vector<NodeOP> NodeOPs;
typedef std::priority_queue<NodeOP, NodeOPs, NodeCompare> NodeQueue;

} // namespace path_finding
} // namespace motif_search

#endif // RNAMAKE_NEW_PATH_FINDING_NODE_H
