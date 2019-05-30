//
// Created by Joseph Yesselman on 3/23/19.
//

#ifndef RNAMAKE_NEW_SOLUTION_TOPOLOGY_H
#define RNAMAKE_NEW_SOLUTION_TOPOLOGY_H

#include <data_structure/graph.h>
#include <motif_data_structure/motif_state_ensemble_graph.h>
#include <resources/motif_state_sqlite_library.h>

namespace motif_search {

class SolutionTopologyTemplate {
public:
    enum class NodeType {
        LIBRARY,
        MOTIF_STATE,
        ENSEMBLE
    };

    class Node {
    public:
        inline
        Node(
                String const & lib_name):
                type_(NodeType::LIBRARY),
                lib_name_(lib_name) {}

        inline
        Node(
                motif::MotifStateOP ms):
                type_(NodeType::MOTIF_STATE),
                motif_state_(ms) {}

        inline
        Node(
                motif::MotifStateEnsembleOP mse):
                type_(NodeType::ENSEMBLE),
                ensemble_(mse) {}

    public:
        inline
        NodeType
        get_type() const { return type_; }

        inline
        String const &
        get_lib_name() const {
            if(type_ != NodeType::LIBRARY) {
                throw std::runtime_error("cannot get_lib_name(), wrong type");
            }
            return lib_name_;
        }

        inline
        motif::MotifStateOP
        get_motif_state() const {
            if(type_ != NodeType::MOTIF_STATE) {
                throw std::runtime_error("cannot get_motif_state(), wrong type");
            }
            return motif_state_;
        }

        inline
        motif::MotifStateEnsembleOP
        get_motif_state_ensemble() const {
            if(type_ != NodeType::ENSEMBLE) {
                throw std::runtime_error("cannot get_motif_state_ensemble(), wrong type");
            }
            return ensemble_;
        }

    private:
        NodeType type_;
        String lib_name_;
        motif::MotifStateOP motif_state_;
        motif::MotifStateEnsembleOP ensemble_;
    };

public:
    SolutionTopologyTemplate():
            g_(data_structure::FixedEdgeDirectedGraph<Node>()) {
        motif_type_ends_ = StringIntMap();
        motif_type_ends_["ideal_helices"] = 2;
        motif_type_ends_["flex_helices"] = 2;
        motif_type_ends_["twoway"]  = 2;
        motif_type_ends_["unique_twoway"] = 2;
    }

    ~SolutionTopologyTemplate() {}

public: // iterators

    typedef typename data_structure::FixedEdgeDirectedGraph<Node>::const_iterator const_iterator;
    typedef typename data_structure::FixedEdgeDirectedGraph<Node>::iterator iterator;

    iterator begin() { return g_.begin(); }
    iterator end()   { return g_.end(); }

    const_iterator begin() const noexcept { return g_.begin(); }
    const_iterator end()   const noexcept { return g_.end(); }

public:
    void
    add_library(
            String const & lib_name) {
        auto num_ends = _get_ends_for_motif_type(lib_name);
        auto n = Node(lib_name);
        g_.add_node(n, num_ends);
        _update_default_transveral();
    }

    void
    add_library(
            String const & lib_name,
            data_structure::NodeIndexandEdge const & parent_nie) {
        auto num_ends = _get_ends_for_motif_type(lib_name);
        auto n = Node(lib_name);
        g_.add_node(n, num_ends, 0, parent_nie);
        _update_default_transveral();
    }

public:
    void
    add_motif_state(
            motif::MotifStateOP ms) {
        auto n = Node(ms);
        g_.add_node(n, ms->end_names().size());
        _update_default_transveral();
    }

    void
    add_motif_state(
            motif::MotifStateOP ms,
            data_structure::NodeIndexandEdge const & parent_nie) {
        auto n = Node(ms);
        g_.add_node(n, ms->end_names().size(), 0, parent_nie);
        _update_default_transveral();
    }

public:
    void
    add_ensemble(
            motif::MotifStateEnsembleOP mse) {
        auto n = Node(mse);
        g_.add_node(n, mse->num_end_states());
        _update_default_transveral();
    }

    void
    add_ensemble(
            motif::MotifStateEnsembleOP mse,
            data_structure::NodeIndexandEdge const & parent_nie) {
        auto n = Node(mse);
        g_.add_node(n, mse->num_end_states(), 0, parent_nie);
        _update_default_transveral();
    }



public:
    inline
    bool
    has_parent(
            Index ni) const {
        return g_.has_parent(ni);
    }

    inline
    Index
    get_parent_index(
            Index ni) const {
        return g_.get_parent_index(ni);
    }

    inline
    Index
    get_parent_end_index(
            Index ni) const {
        return g_.get_parent_end_index(ni);
    }



private:
    void
    _update_default_transveral() {
        auto roots = g_.get_root_indexes();
        if (roots.size() > 0) {
            g_.setup_transversal(roots[0]);
        }
    }

    int
    _get_ends_for_motif_type(
            String const & motif_type) {
        if(motif_type_ends_.find(motif_type) == motif_type_ends_.end()) {
            throw std::runtime_error("motif type: " + motif_type + " not accepted in SolutionTopology");
        }
        return motif_type_ends_[motif_type];
    }

private:
    data_structure::FixedEdgeDirectedGraph<Node> g_;
    StringIntMap motif_type_ends_;

};

typedef std::shared_ptr<SolutionTopologyTemplate> SolutionTopologyTemplateOP;

class SolutionToplogy;
typedef std::shared_ptr<SolutionToplogy> SolutionToplogyOP;

class SolutionToplogyFactory {
public:
    SolutionToplogyFactory()  {}

    ~SolutionToplogyFactory() {}

public:

    SolutionToplogyOP
    generate_toplogy(
            SolutionTopologyTemplate const & sol_template) {

        auto mseg = std::make_shared<motif_data_structure::MotifStateEnsembleOPGraph>();
        for(auto const & n : sol_template) {
            auto mse = motif::MotifStateEnsembleOP();
            if(n->data().get_type() == SolutionTopologyTemplate::NodeType::LIBRARY) {
                auto lib = _get_library(n->data().get_lib_name());
                mse = _parse_library_into_ensemble(lib);
            }
            else if(n->data().get_type() == SolutionTopologyTemplate::NodeType::ENSEMBLE) {
                mse = n->data().get_motif_state_ensemble();
            }
            else if(n->data().get_type() == SolutionTopologyTemplate::NodeType::MOTIF_STATE) {
                mse = std::make_shared<motif::MotifStateEnsemble>(n->data().get_motif_state());
            }

            else {
                throw std::runtime_error("not implemented");
            }

            if(sol_template.has_parent(n->index())) {
                auto parent_nie = data_structure::NodeIndexandEdge{sol_template.get_parent_index(n->index()),
                                                                   sol_template.get_parent_end_index(n->index())};
                mseg->add_ensemble(mse, parent_nie);
            }
            else {
                mseg->add_ensemble(mse);
            }

        }

        return std::make_shared<SolutionToplogy>(mseg);

    }

private:
    resources::MotifStateSqliteLibraryOP
    _get_library(
            String const & lib_name) {
        if(libraries_.find(lib_name) == libraries_.end()) {
            libraries_[lib_name] = std::make_shared<resources::MotifStateSqliteLibrary>(lib_name);
            libraries_[lib_name]->load_all();
        }
        return libraries_[lib_name];
    }

    motif::MotifStateEnsembleOP
    _parse_library_into_ensemble(
            resources::MotifStateSqliteLibraryOP library) {
        auto motif_states = motif::MotifStateOPs();
        auto energies = Floats();

        for(auto const & ms : *library) {
            motif_states.push_back(ms);
            energies.push_back(1);
        }

        return std::make_shared<motif::MotifStateEnsemble>(motif_states, energies);
    }



private:
    base::Options options_;
    std::map<String, resources::MotifStateSqliteLibraryOP> libraries_;
    std::map<String, motif::MotifStateEnsembleOP> ensembles_;

};

class SolutionToplogy {
public:
    SolutionToplogy(
            motif_data_structure::MotifStateEnsembleOPGraphOP mseg) {
        mseg_ = mseg;
        rng_ = util::RandomNumberGenerator();
        solution_nie_ = mseg_->get_leafs();
    }

public:

    typedef typename motif_data_structure::MotifStateEnsembleOPGraph::const_iterator const_iterator;
    typedef typename motif_data_structure::MotifStateEnsembleOPGraph::iterator iterator;

    iterator begin() { return mseg_->begin(); }
    iterator end()   { return mseg_->end(); }

    const_iterator begin() const noexcept { return mseg_->begin(); }
    const_iterator end()   const noexcept { return mseg_->end(); }

public:
    motif_data_structure::MotifStateGraphOP
    initialize_solution(
            structure::BasepairStateOP bp_state) {

        auto ms = std::make_shared<motif::MotifState>(
                "start", Strings{"start", "start"}, Strings{"", ""},
                structure::BasepairStateOPs {bp_state, bp_state},
                math::Points(), 0, 0, 0);

        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg->set_option_value("sterics", false);
        msg->add_state(ms);

        for (auto const & n : *mseg_) {
            auto ms = get_motif_state(n->index());
            if (mseg_->has_parent(n->index())) {
                msg->add_state(ms, mseg_->get_parent_index(n->index()) + 1, mseg_->get_parent_end_index(n->index()));
            } else {
                msg->add_state(ms);
            }

        }
        return msg;
    }

    motif::MotifStateOP
    get_motif_state(
            Index pos) {
        max_member_ = mseg_->get_ensemble(pos)->size();
        mem_pos_ = rng_.randrange(max_member_);
        return mseg_->get_ensemble(pos)->get_member(mem_pos_)->motif_state;
    }

    std::vector<data_structure::NodeIndexandEdge> const &
    get_solution_nie() {
        return solution_nie_;
    }

    inline
    size_t
    size() {
        return mseg_->size();
    }

    inline
    size_t
    get_ensemble_size(
          Index pos) {
        return mseg_->get_ensemble(pos)->size();
    }

private:
    motif_data_structure::MotifStateEnsembleOPGraphOP mseg_;
    util::RandomNumberGenerator rng_;
    int max_member_, mem_pos_;
    std::vector<data_structure::NodeIndexandEdge> solution_nie_;
};

}

#endif //RNAMAKE_NEW_SOLUTION_TOPOLOGY_H


























