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
            float score):
            state_(ms),
            parent_(nullptr),
            score_(score),
            ss_score_(){}

    inline
    Node(
            motif::MotifStateOP ms,
            NodeOP parent,
            float score):
            state_(ms),
            parent_(parent),
            score_(score){}


    ~Node() {}


public: // getters

    int
    level() { return level_; }

    inline
    int
    size() { return size_; }

    inline
    int
    ntype() { return ntype_; }

    inline
    float
    ss_score() { return ss_score_; }

    inline
    float
    score() const { return score_; }

    inline
    int
    parent_end_index() { return parent_end_index_; }

    inline
    motif::MotifStateOP
    state() const { return state_;}

private:
    NodeOP parent_;
    motif::MotifStateOP state_;
    int parent_end_index_, level_, size_, ntype_;
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
