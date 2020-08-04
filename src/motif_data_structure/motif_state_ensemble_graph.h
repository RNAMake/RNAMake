//
// Created by Joseph Yesselman on 4/25/18.
//

#ifndef TEST_MOTIF_STATE_ENSEMBLE_GRAPH_H
#define TEST_MOTIF_STATE_ENSEMBLE_GRAPH_H

#include <data_structure/graph.h>
#include "data_structure/graph/graph.h"
#include "motif/motif_state_ensemble.h"
#include "motif_data_structure/motif_graph.h"
#include "motif_data_structure/motif_state_graph.hpp"


namespace motif_data_structure {

class MotifStateEnsembleGraphException : public std::runtime_error {
public:
    MotifStateEnsembleGraphException(
            String const & message) :
            std::runtime_error(message) {}
};

class MotifStateEnsembleGraph {
public:

    MotifStateEnsembleGraph():
            graph_(data_structure::FixedEdgeDirectedGraph<motif::MotifStateEnsemble>()) {
    }

    MotifStateEnsembleGraph(
            MotifStateEnsembleGraph const & mseg) {
        graph_ = mseg.graph_;
        _update_default_transveral();
    }

public: // iterator

    typedef typename data_structure::FixedEdgeDirectedGraph<motif::MotifStateEnsemble>::const_iterator const_iterator;
    typedef typename data_structure::FixedEdgeDirectedGraph<motif::MotifStateEnsemble>::iterator iterator;

    iterator begin() { return graph_.begin(); }
    iterator end()   { return graph_.end(); }

    const_iterator begin() const noexcept { return graph_.begin(); }
    const_iterator end()   const noexcept { return graph_.end(); }

public:

    size_t
    size() { return graph_.get_num_nodes(); }

public: // add functions

    int
    add_ensemble(
            motif::MotifStateEnsemble const & ensemble) {
        auto ni = graph_.add_node(ensemble, ensemble.num_end_states());
        _update_default_transveral();
        return ni;
    }

    int
    add_ensemble(
            motif::MotifStateEnsemble const & ensemble,
            data_structure::NodeIndexandEdge const & parent_nie) {
        auto ni = graph_.add_node(ensemble, ensemble.num_end_states(), 0, parent_nie);
        _update_default_transveral();
        return ni;
    }

public: // get

    inline
    motif::MotifStateEnsemble const &
    get_ensemble(
            Index ni) {
        return graph_.get_node_data(ni);
    }

public:
    inline
    bool
    has_parent(
            Index ni) const {
        return graph_.has_parent(ni);
    }

    inline
    Index
    get_parent_index(
            Index ni) const {
        return graph_.get_parent_index(ni);
    }

    inline
    Index
    get_parent_end_index(
            Index ni) const {
        return graph_.get_parent_end_index(ni);
    }

public: // getters for connections

    inline
    std::vector<data_structure::Edge const *> const &
    get_ensemble_connections(
            Index ni) const {
        return graph_.get_node_edges(ni);
    }

    inline
    bool
    are_ensembles_connected(
            Index n1,
            Index n2) {
        return graph_.edge_between_nodes(n1, n2);
    }

private:
    void
    _update_default_transveral() {
        auto roots = graph_.get_root_indexes();
        if (roots.size() > 0) {
            graph_.setup_transversal(roots[0]);
        }
    }

private:
    data_structure::FixedEdgeDirectedGraph<motif::MotifStateEnsemble> graph_;

};

class MotifStateEnsembleOPGraph {
public:

    MotifStateEnsembleOPGraph():
            graph_(data_structure::FixedEdgeDirectedGraph<motif::MotifStateEnsembleOP>()) {
    }

public: // iterator

    typedef typename data_structure::FixedEdgeDirectedGraph<motif::MotifStateEnsembleOP>::const_iterator const_iterator;
    typedef typename data_structure::FixedEdgeDirectedGraph<motif::MotifStateEnsembleOP>::iterator iterator;

    iterator begin() { return graph_.begin(); }
    iterator end()   { return graph_.end(); }

    const_iterator begin() const noexcept { return graph_.begin(); }
    const_iterator end()   const noexcept { return graph_.end(); }

public:

    size_t
    size() { return graph_.get_num_nodes(); }

public: // add functions

    int
    add_ensemble(
            motif::MotifStateEnsembleOP ensemble) {
        auto ni = graph_.add_node(ensemble, ensemble->num_end_states());
        _update_default_transveral();
        return ni;
    }

    int
    add_ensemble(
            motif::MotifStateEnsembleOP ensemble,
            data_structure::NodeIndexandEdge const & parent_nie) {
        auto ni = graph_.add_node(ensemble, ensemble->num_end_states(), 0, parent_nie);
        _update_default_transveral();
        return ni;
    }

public: // get functions

    inline
    motif::MotifStateEnsembleOP
    get_ensemble(
            Index ni) {
        return graph_.get_node_data(ni);
    }

public:
    inline
    bool
    has_parent(
            Index ni) const {
        return graph_.has_parent(ni);
    }

    inline
    Index
    get_parent_index(
            Index ni) const {
        return graph_.get_parent_index(ni);
    }

    inline
    Index
    get_parent_end_index(
            Index ni) const {
        return graph_.get_parent_end_index(ni);
    }

    inline
    std::vector<data_structure::NodeIndexandEdge>
    get_leafs() {
        auto leafs = std::vector<data_structure::NodeIndexandEdge>();
        for(auto const & n : graph_) {
            auto num_edges = n->data()->num_end_states();
            for(int i = 1 ; i < num_edges; i++) {
                if(graph_.edge_index_empty(n->index(), i)) {
                    leafs.push_back(data_structure::NodeIndexandEdge{n->index(), i});
                }
            }
        }
        return leafs;
    }


private:
    void
    _update_default_transveral() {
        auto roots = graph_.get_root_indexes();
        if (roots.size() > 0) {
            graph_.setup_transversal(roots[0]);
        }
    }

private:
    data_structure::FixedEdgeDirectedGraph<motif::MotifStateEnsembleOP> graph_;

};


typedef std::shared_ptr<MotifStateEnsembleGraph> MotifStateEnsembleGraphOP;
typedef std::shared_ptr<MotifStateEnsembleOPGraph> MotifStateEnsembleOPGraphOP;

}
#endif //TEST_MOTIF_STATE_ENSEMBLE_GRAPH_H
