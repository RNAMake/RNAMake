//
//  motif_graph.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_graph__
#define __RNAMake__motif_graph__

#include <stdio.h>
#include <map>

//RNAMake Headers
#include "base/option.h"
#include "data_structure/graph/graph.h"
#include "motif/motif.h"
#include "motif_data_structure/motif_tree.h"
#include "motif_data_structure/motif_merger.h"

namespace motif_data_structure {

enum MotifGraphStringType {
    OLD,
    MG,
    TOP
};

class MotifGraphException : public std::runtime_error {
public:
    MotifGraphException(
            String const & message) :
            std::runtime_error(message) {}
};

class MotifGraph {

public: //construction

    MotifGraph();

    MotifGraph(
            String const &);

    MotifGraph(
            String const &,
            MotifGraphStringType const &);

    MotifGraph(
            MotifGraph const &);

    ~MotifGraph() {}

public: //setup helpers

    void
    update_indexes(std::map<int, int> const &);

public:

    typedef typename data_structure::graph::GraphNodeOPs<motif::MotifOP>::const_iterator const_iterator;
    typedef typename data_structure::graph::GraphNodeOPs<motif::MotifOP>::iterator iterator;

    iterator begin() {
        _update_align_list();
        return align_list_.begin();
    }

    iterator end() {
        return align_list_.end();
    }


private://add function helpers

    data_structure::graph::GraphNodeOP<motif::MotifOP>
    _get_parent(
            String const &,
            int);

    int
    _add_motif_to_graph(
            motif::MotifOP &,
            data_structure::graph::GraphNodeOP<motif::MotifOP> const &,
            int);

    Ints
    _get_available_parent_end_pos(
            data_structure::graph::GraphNodeOP<motif::MotifOP> const &,
            int);

    void
    _add_motif_tree(
            MotifTreeOP const &,
            int,
            int);

    int
    _steric_clash(
            motif::MotifOP const &);


    int
    _get_connection_end(
            data_structure::graph::GraphNodeOP<motif::MotifOP> const &,
            String const &);

    void
    _align_motifs_all_motifs();

public: //add functions

    int
    add_motif(
            motif::MotifOP const & m,
            int parent_index = -1,
            int parent_end_index = -1,
            int orphan = 0);

    int
    add_motif(
            motif::MotifOP const &,
            int,
            String const &);

    void
    add_motif_tree(
            MotifTreeOP const & mt,
            int parent_index = -1,
            int parent_end_index = -1);

    void
    add_motif_tree(
            MotifTreeOP const &,
            int,
            String const &);

    void
    add_connection(
            int,
            int,
            String const &,
            String const &);

public: //remove functions

    void
    remove_motif(int);

    void
    remove_level(int level);


public: //graph wrappers
    inline
    size_t
    size() { return graph_.size(); }

    inline
    void
    increase_level() { return graph_.increase_level(); }

    inline
    data_structure::graph::GraphNodeOP<motif::MotifOP> const &
    last_node() { return graph_.last_node(); }

    inline
    data_structure::graph::GraphNodeOP<motif::MotifOP>
    oldest_node() { return graph_.oldest_node(); }

    inline
    data_structure::graph::GraphNodeOP<motif::MotifOP> const &
    get_node(int i) const { return graph_.get_node(i); }

    inline
    data_structure::graph::GraphNodeOP<motif::MotifOP> const
    get_node(util::Uuid const & uuid) const {
        for (auto const & n : graph_) {
            if (n->data()->id() == uuid) {
                return n;
            }
        }
        throw MotifGraphException("cannot get node with uuid no motif has it in this graph");
    }

    inline
    data_structure::graph::GraphNodeOP<motif::MotifOP> const
    get_node(String const & m_name) const {
        auto node = data_structure::graph::GraphNodeOP<motif::MotifOP>(nullptr);
        for (auto const & n : graph_) {
            if (n->data()->name() == m_name) {
                if (node != nullptr) {
                    throw MotifGraphException(
                            "cannot get node with name: " + m_name + " there is more then one motif "
                                    "with this name");
                }

                node = n;
            }
        }

        if (node == nullptr) {
            throw MotifGraphException(
                    "cannot get node with name: " + m_name + " there is no motif in the tree with "
                            "this name");
        }

        return node;
    }

    inline
    void
    set_index(int index) { graph_.index(index); }


public: //designing functions
    void
    replace_ideal_helices();

    void
    replace_helical_sequence(
            secondary_structure::PoseOP const &);

    inline
    void
    replace_helical_sequence(
            String const & seq) {
        auto dss = designable_secondary_structure();
        dss->replace_sequence(seq);
        replace_helical_sequence(dss);
    }

    secondary_structure::PoseOP
    designable_secondary_structure() {
        _update_merger();
        auto ss = merger_->secondary_structure();
        auto ss_r = secondary_structure::ResidueOP(nullptr);

        for (auto const & n : graph_) {
            if (n->data()->name() != "HELIX.IDEAL") { continue; }

            for (auto const & r : n->data()->residues()) {
                ss_r = ss->get_residue(r->uuid());
                if (ss_r != nullptr) {
                    ss_r->name("N");
                }
            }
        }

        return ss;
    }

    inline
    String
    designable_sequence() {
        return designable_secondary_structure()->sequence();
    }


public: // outputing functions
    void
    write_pdbs(String const & fname = "nodes");

    String
    topology_to_str();

    String
    to_str();

public: // misc functions

    void
    _update_align_list();

    void
    _update_merger();

public: // getters

    inline
    structure::Beads
    beads() {
        structure::Beads beads;
        for (auto const & n : graph_.nodes()) {
            std::copy(n->data()->beads().begin(),
                      n->data()->beads().end(),
                      std::inserter(beads, beads.end()));
        }
        return beads;
    }

    structure::BasepairOP const &
    get_available_end(int);

    structure::BasepairOP const &
    get_available_end(
            int,
            String const &);

    structure::BasepairOP
    get_available_end(
            String const &,
            String const &);

    data_structure::graph::GraphNodeOPs<motif::MotifOP> const
    unaligned_nodes() const;

    // this is bad
    std::map<int, int>
    aligned() { return aligned_; }

    data_structure::graph::GraphConnectionOPs<motif::MotifOP> const &
    connections() { return graph_.connections(); }


public: //Motif Merger Wrappers

    inline
    structure::RNAStructureOP const &
    get_structure() {
        try {
            _update_merger();
            return merger_->get_structure();
        }
        catch (MotifMergerException) {
            throw MotifGraphException(
                    "cannot produce merged structure it is likely you have created a ring with no start"
                            "call write_pdbs() to see what the topology would look like");
        }
    }

    secondary_structure::PoseOP
    secondary_structure() {

        _update_merger();
        return merger_->secondary_structure();

        try {
            _update_merger();
            return merger_->secondary_structure();
        }
        catch (MotifMergerException) {
            throw MotifGraphException(
                    "cannot produce merged secondary structure it is likely you have created a ring "
                            "with no start, call write_pdbs() to see what the topology would look like");
        }

    }

    inline
    void
    to_pdb(
            String const fname = "test.pdb",
            int renumber = -1,
            int close_chains = 0,
            int conect_statement = 0) {

        try {
            _update_merger();
            return merger_->to_pdb(fname, renumber, close_chains, conect_statement);
        }
        catch (MotifMergerException) {
            throw MotifGraphException(
                    "cannot produce merged structure for a pdb it is likely you have created a ring "
                            "with no start, call write_pdbs() to see what the topology would look like");
        }

    }

    inline
    String
    sequence() {
        return secondary_structure()->sequence();
    }

    inline
    String
    dot_bracket() {
        return secondary_structure()->dot_bracket();
    }


public: //Options Wrappers

    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }

    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }

    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }

    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }


    template<typename T>
    void
    set_option_value(
            String const & name,
            T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }


private:
    void
    setup_options();

    void
    update_var_options();

private:
    void
    _setup_from_top_str(String const &);

    void
    _setup_from_str(String const &);

private:
    data_structure::graph::GraphStatic<motif::MotifOP> graph_;
    data_structure::graph::GraphNodeOPs<motif::MotifOP> align_list_;
    MotifMergerOP merger_;
    base::Options options_;
    std::map<int, int> aligned_;
    int update_merger_;
    int update_align_list_;
    //options
    float clash_radius_;
    bool sterics_;

};

typedef std::shared_ptr<MotifGraph> MotifGraphOP;

}

#endif /* defined(__RNAMake__motif_graph__) */
