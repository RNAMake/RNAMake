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

    inline
    Node(
            motif::MotifStateOP ms,
            NodeOP parent,
            float score,
            int level,
            int parent_end_index,
            int node_type):
            state_(ms),
            parent_(parent),
            score_(score),
            level_(level),
            parent_end_index_(parent_end_index),
            node_type_(node_type),
            size_(ms->size()),
            ss_score_(ms->score()) {
        if(parent_ != nullptr) {
            ss_score_ += parent_->ss_score();
            size_ += parent->size();
        }
    }

    ~Node() {}


public: // getters

    int
    level() const { return level_; }

    inline
    int
    size() const { return size_; }

    inline
    int
    node_type() const { return node_type_; }

    inline
    float
    ss_score() const { return ss_score_; }

    inline
    float
    score() const { return score_; }

    inline
    int
    parent_end_index() const { return parent_end_index_; }

    inline
    motif::MotifStateOP
    state() const { return state_;}

    inline
    NodeOP
    parent() const { return parent_; }

private:
    NodeOP parent_;
    motif::MotifStateOP state_;
    int parent_end_index_, level_, size_, node_type_;
    float ss_score_, score_;
};

struct NodeCompare {
    bool
    operator()(
            NodeOP node1,
            NodeOP node2) {

        if (node1->score() > node2->score()) { return true; }
        else { return false; }
    }
};

typedef std::vector<NodeOP> NodeOPs;
typedef std::priority_queue<NodeOP, NodeOPs, NodeCompare> NodeQueue;

}
}


#endif //RNAMAKE_NEW_PATH_FINDING_NODE_H
