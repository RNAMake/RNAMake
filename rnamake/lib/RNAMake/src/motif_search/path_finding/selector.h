//
// Created by Joseph Yesselman on 3/18/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_SELECTOR_H
#define RNAMAKE_NEW_PATH_FINDING_SELECTOR_H

#include <data_structure/graph/graph.h>
#include <motif/motif_state.h>
#include <motif/motif_state_ensemble.h>
#include <resources/motif_state_sqlite_library.h>

namespace motif_search {
namespace path_finding {

class SelectorException : public std::runtime_error {
public:
    SelectorException(
            String const & message) :
            std::runtime_error(message) {}
};

struct SelectorNodeData {
    inline
    SelectorNodeData(
            String const & n_name,
            motif::MotifStateOPs const & n_motif_states,
            int n_type):
            name(n_name),
            motif_states(n_motif_states),
            type(n_type) {}

    String name;
    motif::MotifStateOPs motif_states;
    int type;
};

typedef std::shared_ptr<SelectorNodeData> SelectorNodeDataOP;

class Selector {
public:
    Selector():
            pos_(0),
            max_(0),
            parent_type_(0) {}

    virtual
    ~Selector() {}

    virtual
    Selector *
    clone() const {
        return new Selector(*this);
    }

public:

    virtual
    void
    add (
            String const & lib_name) {
        auto type = graph_.size();
        auto motif_states = motif::MotifStateOPs();
        auto ms_lib = resources::MotifStateSqliteLibrary(lib_name);
        ms_lib.load_all();
        for (auto const & ms : ms_lib) { motif_states.push_back(ms); }
        auto d = std::make_shared<SelectorNodeData>(lib_name, motif_states, type);
        _add(d);
    }

    virtual
    void
    add (
            String const & lib_name,
            motif::MotifStateEnsembleOP mse) {
        auto type = graph_.size();
        auto motif_states = motif::MotifStateOPs();
        for (auto const & mem : mse->members()) { motif_states.push_back(mem->motif_state); }
        auto d = std::make_shared<SelectorNodeData>(lib_name, motif_states, type);
        _add(d);
    }

    virtual
    void
    add (
            motif::MotifOP motif) {
        auto type = graph_.size();
        auto motif_states = motif::MotifStateOPs{motif->get_state()};
        auto d = std::make_shared<SelectorNodeData>(motif->name(), motif_states, type);
        _add(d);
    }

public:

    void
    connect(
            String const &,
            String const &);

public:
    void
    start(
            int parent_type) const {
        parent_type_ = parent_type;
        if(parent_type == -1) {
            max_ = 0;
        }
        else {
            max_ = graph_.get_node(parent_type)->connections().size() - 1;
        }
        pos_ = -1;

    }

    SelectorNodeDataOP
    next() const {
        if(finished()) {
            throw SelectorException("cannot enumerate anymore selector is done enumerating");
        }
        pos_ += 1;

        if(parent_type_ == -1) {
            return graph_.get_node(0)->data();
        }

        else {
            auto c = graph_.get_node(parent_type_)->connections()[pos_];
            std::cout << c->end_index_1() << " " << c->end_index_2() << " " << c->node_1() << " " << c->node_2() << std::endl;
            return c->partner(parent_type_)->data();
        }

    }

    bool
    finished() const {
        return pos_ == max_;
    }

    size_t
    size() const { return graph_.size(); }


protected:
    void
    _add(
            SelectorNodeDataOP d) {
        for(auto const & n : graph_) {
            if(d->name == n->data()->name) {
                throw SelectorException("cannot have two nodes with the same name");
            }
        }
        graph_.add_data(d, -1, 1);

    }

protected:
    data_structure::graph::GraphDynamic<SelectorNodeDataOP> graph_;
    mutable int pos_, max_, parent_type_;

};

class RoundRobinSelector : public Selector {
public:
    RoundRobinSelector(): Selector() {}

    ~RoundRobinSelector() {}

    Selector *
    clone() const {
        return new RoundRobinSelector(*this);
    }

public:
    virtual
    void
    add (
            String const & lib_name) {
        Selector::add(lib_name);
        _add_new_connections();
    }

    virtual
    void
    add (
            String const & lib_name,
            motif::MotifStateEnsembleOP mse) {
        Selector::add(lib_name, mse);
        _add_new_connections();
    }

    virtual
    void
    add (
            motif::MotifOP motif) {
        Selector::add(motif);
        _add_new_connections();
    }

private:
    void
    _add_new_connections() {
        auto i = (int) graph_.size() - 1;
        for (int j = 0; j < graph_.size(); j++) {
            if(j == i) { continue; }
            graph_.connect(j, i);
        }
    }
};

typedef std::shared_ptr<Selector> SelectorOP;

SelectorOP
default_selector();

}
}


#endif //RNAMAKE_NEW_PATH_FINDING_SELECTOR_H
