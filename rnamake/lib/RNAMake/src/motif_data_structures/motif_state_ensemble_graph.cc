//
// Created by Joseph Yesselman on 4/25/18.
//

#include "motif_state_ensemble_graph.h"


MotifStateEnsembleGraph::MotifStateEnsembleGraph():
        graph_(GraphStatic<MotifStateEnsembleOP>()),
        update_align_list_(1),
        align_list_(GraphNodeOPs<MotifStateEnsembleOP>()),
        aligned_(std::map<int, int>()) { }


MotifStateEnsembleGraph::MotifStateEnsembleGraph(
        MotifGraphOP const & mg):
        MotifStateEnsembleGraph() {

}


//add function helpers /////////////////////////////////////////////////////////////////////////////


GraphNodeOP<MotifStateEnsembleOP>
MotifStateEnsembleGraph::_get_parent(
        int parent_index) {

    auto parent = graph_.last_node();

    //catch non existant parent
    try {
        if(parent_index != -1) { parent = graph_.get_node(parent_index); }
    }
    catch(GraphException const & e) {
        throw MotifStateEnsembleGraphException (
                "could not add state ensemble with parent index: " +
                std::to_string(parent_index) + "there is no node with that index");
    }

    return parent;

}

Ints
MotifStateEnsembleGraph::_get_available_parent_end_pos(
        GraphNodeOP<MotifStateEnsembleOP> const & parent,
        int parent_end_index) {

    auto avail_pos = Ints();

    if(parent_end_index != -1) {
        int avail = parent->available_pos(parent_end_index);
        if(!avail) {
            throw MotifStateEnsembleGraphException(
                    "cannot add state ensemble to graph as the end with index: " +
                    std::to_string(parent_end_index) +" you are trying to "
                            "add it to is already filled or does not exist");
        }

        if(parent_end_index == parent->data()->block_end_add()) {
            throw MotifStateEnsembleGraphException(
                    "cannot add motif: to tree as the parent_end_index"
                            " supplied is blocked see class Motif");
        }

        avail_pos.push_back(parent_end_index);
    }

    else {
        auto avail_pos_temp = parent->available_children_pos();
        for(auto const & p : avail_pos_temp) {
            if(p == parent->data()->block_end_add()) { continue; }
            avail_pos.push_back(p);
        }
    }

    return avail_pos;

}


//add function  ////////////////////////////////////////////////////////////////////////////////////


int
MotifStateEnsembleGraph::add_ensemble(
        MotifStateEnsembleOP const & ensemble,
        int parent_index,
        int parent_end_index) {

    auto parent = _get_parent(parent_index);
    if(parent == nullptr) {
        auto pos = graph_.add_data(std::make_shared<MotifStateEnsemble>(*ensemble), -1, -1, -1,
                                   (int)ensemble->num_end_states());
        return pos;
    }

    auto avail_pos = _get_available_parent_end_pos(parent, parent_end_index);
    if(avail_pos.size() == 0) {
        throw MotifStateEnsembleGraphException("no available positions to add too");
    }
    if(avail_pos.size() > 1) {
        throw MotifStateEnsembleGraphException("more than one positions to add ensemble too");
    }

    auto pos = graph_.add_data(std::make_shared<MotifStateEnsemble>(*ensemble), parent->index(), avail_pos[0], 0,
                               (int)ensemble->num_end_states());
    if(pos != -1) {
        update_align_list_ = 1;
        aligned_[pos] = 1;
    }
    return pos;
}


void
MotifStateEnsembleGraph::add_connection(
        int ni,
        int nj,
        int ei,
        int ej) {

    graph_.connect(ni, nj, ei, ej);
}

