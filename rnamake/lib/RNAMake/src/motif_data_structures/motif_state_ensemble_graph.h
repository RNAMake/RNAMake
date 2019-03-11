//
// Created by Joseph Yesselman on 4/25/18.
//

#ifndef TEST_MOTIF_STATE_ENSEMBLE_GRAPH_H
#define TEST_MOTIF_STATE_ENSEMBLE_GRAPH_H

#include "data_structure/graph/graph.h"
#include "motif/motif_state_ensemble.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_data_structures/motif_state_graph.hpp"


class MotifStateEnsembleGraphException : public std::runtime_error {
public:
    MotifStateEnsembleGraphException(
            String const & message):
            std::runtime_error(message)
    {}
};

class MotifStateEnsembleGraph {
public:

    MotifStateEnsembleGraph();

    MotifStateEnsembleGraph(MotifGraphOP const &);

    MotifStateEnsembleGraph(MotifStateGraphOP const &);

    // copy constuctor
    MotifStateEnsembleGraph(MotifStateEnsembleGraph const &);

public: // iterator

    typedef typename data_structure::graph::GraphStatic<MotifStateEnsembleOP>::iterator iterator;
    typedef typename data_structure::graph::GraphStatic<MotifStateEnsembleOP>::const_iterator const_iterator;

    iterator begin() { return graph_.begin(); }
    iterator end()   { return graph_.end(); }

    const_iterator begin() const { return graph_.begin(); }
    const_iterator end()   const { return graph_.end(); }

public:

    size_t
    size() { return graph_.size(); }

private: //add function helpers

    data_structure::graph::GraphNodeOP<MotifStateEnsembleOP>
    _get_parent(
            int);


    Ints
    _get_available_parent_end_pos(
            data_structure::graph::GraphNodeOP<MotifStateEnsembleOP> const &,
            int);


public: // add functions

    int
    add_ensemble(
            MotifStateEnsembleOP const & ensemble,
            int parent_index = -1,
            int parent_end_index = -1);

    void
    add_connection(
            int,
            int,
            int,
            int);


private:
    data_structure::graph::GraphStatic<MotifStateEnsembleOP> graph_;
    data_structure::graph::GraphNodeOPs<MotifStateEnsembleOP> align_list_;
    std::map<int, int> aligned_;
    int update_align_list_;

};


#endif //TEST_MOTIF_STATE_ENSEMBLE_GRAPH_H
