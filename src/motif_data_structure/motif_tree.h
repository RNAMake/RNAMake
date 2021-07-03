//
//  motif_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree__
#define __RNAMake__motif_tree__

#include <stdio.h>
#include <iomanip>      // std::setw
#include "base/option.h"
#include "data_structure/tree/tree.h"
#include "motif/motif.h"
#include "motif_data_structure/motif_merger.h"
#include "motif_data_structure/motif_connection.h"

namespace motif_data_structure {

enum MotifTreeStringType {
    MT_STR,
    TOP_STR
};

class MotifTreeException : public std::runtime_error {
public:
    MotifTreeException(
            String const & message) :
            std::runtime_error(message) {}
};

class MotifTree {
private:

    class MotifTreePrinter {
    public:
        inline
        MotifTreePrinter(
                MotifTree const & mt) :
                branch_length_(25),
                levels_(std::map<int, int>()),
                node_pos_(std::map<int, int>()),
                nodes_per_level_(std::map<int, int>()),
                start_pos_(100) {
            _setup_node_positions(mt);
        }

        ~MotifTreePrinter() {}

    public:
        String
        print_tree(
                MotifTree const & mt) {
            auto nodes_per_level = std::map<int, data_structure::tree::TreeNodeOPs<motif::MotifOP>>();
            for (auto const & n : mt) {
                nodes_per_level[levels_[n->index()]].push_back(n);
            }

            int found = 1, level = 1;
            auto s = String("\n");

            while (found) {
                if (nodes_per_level.find(level) == nodes_per_level.end()) { break; }
                auto node_level = nodes_per_level[level];
                s += _print_level(node_level);
                level += 1;
            }
            return s;
        }

    private:
        String
        _print_level(data_structure::tree::TreeNodeOPs<motif::MotifOP> const & nodes) {
            auto s = String("");
            auto strings = std::vector<Strings>();

            //sort nodes by pos, using a lambda statement and std::pair
            //probably can be cleaned up
            auto pairs = std::vector<std::pair<data_structure::tree::TreeNodeOP<motif::MotifOP>, int>>();
            for (auto const & n : nodes) {
                auto p = std::pair<data_structure::tree::TreeNodeOP<motif::MotifOP>, int>(n, node_pos_[n->index()]);
                pairs.push_back(p);
            }

            std::sort(pairs.begin(), pairs.end(),
                      [](std::pair<data_structure::tree::TreeNodeOP<motif::MotifOP>, int> const & p1,
                              std::pair<data_structure::tree::TreeNodeOP<motif::MotifOP>, int> const & p2) {
                          return p1.second < p2.second;
                      });


            auto sorted_nodes = data_structure::tree::TreeNodeOPs<motif::MotifOP>();
            for (auto const & p : pairs) { sorted_nodes.push_back(p.first); }

            for (auto const & n : sorted_nodes) {
                auto strs = Strings();
                if (n->parent() != nullptr) {
                    auto parent_end_index = n->parent_end_index();
                    auto parent_end_name = n->parent()->data()->ends()[parent_end_index]->name();
                    strs.push_back("|");
                    strs.push_back("E" + std::to_string(parent_end_index) + " - " + parent_end_name);
                    strs.push_back("|");
                }

                strs.push_back("N" + std::to_string(n->index()) + " - " + n->data()->name());
                strs.push_back("|  - " + n->data()->ends()[0]->name());
                strings.push_back(strs);
            }

            auto transposed_strings = std::vector<Strings>();
            for (int i = 0; i < strings[0].size(); i++) {
                auto transposed = Strings();
                for (auto const & strs: strings) {
                    transposed.push_back(strs[i]);
                }
                transposed_strings.push_back(transposed);
            }

            for (auto const & strs : transposed_strings) {
                int current_pos = 0;
                int j = 0;
                for (auto const & n : sorted_nodes) {
                    int pos = node_pos_[n->index()];
                    int diff = pos - current_pos;
                    for (int i = 0; i < diff; i++) { s += " "; }
                    s += strs[j];
                    current_pos = pos + strs[j].length();
                    j += 1;
                }
                s += "\n";
            }

            int current_pos = 0;
            int hit = 0;
            for (auto const & n : sorted_nodes) {
                auto children = data_structure::tree::TreeNodeOPs<motif::MotifOP>();
                for (auto const & c : n->children()) {
                    if (c != nullptr) { children.push_back(c); }
                }

                if (children.size() > 1) {
                    hit = 1;
                    auto min = node_pos_[children[0]->index()];
                    auto max = node_pos_[children.back()->index()];
                    auto diff = min + 1 - current_pos;
                    for (int i = 0; i < diff; i++) { s += " "; }
                    for (int i = 0; i < max - min - 1; i++) { s += "_"; }
                    current_pos = max;
                }
            }

            if (hit) { s += "\n"; }

            return s;
        }

        void
        _setup_node_positions(
                MotifTree const & mt) {

            _assign_node_levels(mt);

            for (auto const & n : mt) {
                auto n_level = levels_[n->index()];
                if (nodes_per_level_.find(n_level) == nodes_per_level_.end()) {
                    nodes_per_level_[n_level] = 0;
                }
                nodes_per_level_[n_level] += 1;
            }

            int i = -1;
            for (auto const & n : mt) {
                i++;
                if (i == 0) {
                    node_pos_[n->index()] = start_pos_;
                }

                auto children = data_structure::tree::TreeNodeOPs<motif::MotifOP>();
                for (auto const & c : n->children()) {
                    if (c != nullptr) { children.push_back(c); }
                }

                if (children.size() == 1) {
                    node_pos_[children[0]->index()] = node_pos_[n->index()];
                } else if (children.size() == 2) {
                    auto level = levels_[n->index()];
                    auto nodes_per_level = nodes_per_level_[level + 1];
                    auto extra = nodes_per_level - 2;
                    auto parent_pos = node_pos_[n->index()];
                    if (extra == 0) {
                        node_pos_[children[0]->index()] = parent_pos - branch_length_;
                        node_pos_[children[1]->index()] = parent_pos + branch_length_;
                    } else {
                        node_pos_[children[0]->index()] = parent_pos - branch_length_ / extra;
                        node_pos_[children[1]->index()] = parent_pos + branch_length_ / extra;
                    }
                } else if (children.size() > 2) {
                    throw MotifTreeException(
                            "Greater then two children is not supported for pretty_printing");
                }

            }
        }

        void
        _assign_node_levels(
                MotifTree const & mt) {

            for (auto const & n : mt) {
                if (levels_.size() == 0) {
                    levels_[n->index()] = 1;
                } else {
                    auto parent_level = levels_[n->parent_index()];
                    levels_[n->index()] = parent_level + 1;
                }
            }
        }


    private:
        std::map<int, int> levels_;
        std::map<int, int> node_pos_;
        std::map<int, int> nodes_per_level_;
        int branch_length_;
        int start_pos_;


    };

public: // constructors 

    MotifTree();

    MotifTree(
            String const &);

    MotifTree(
            String const &,
            MotifTreeStringType);

    MotifTree(
            MotifTree const &);

    ~MotifTree() {}

private: // setup functions

    void
    _setup_from_str(
            String const &);

public: //iterators

    typedef typename data_structure::tree::TreeStatic<motif::MotifOP>::iterator iterator;
    typedef typename data_structure::tree::TreeStatic<motif::MotifOP>::const_iterator const_iterator;

    iterator begin() { return tree_.begin(); }

    iterator end() { return tree_.end(); }

    const_iterator begin() const { return tree_.begin(); }

    const_iterator end() const { return tree_.end(); }

private: //adding functions helpers

    data_structure::tree::TreeNodeOP<motif::MotifOP>
    _get_parent(int);

    Ints
    _get_available_parent_end_pos(
            data_structure::tree::TreeNodeOP<motif::MotifOP> const &,
            int);

    int
    _steric_clash(
            motif::MotifOP const &);

    int
    _get_connection_end(
            data_structure::tree::TreeNodeOP<motif::MotifOP> const &,
            String const &);


public: //adding functions


    int
    add_motif(
            motif::MotifOP const & m,
            int parent_index = -1,
            int parent_end_index = -1);

    int
    add_motif(
            motif::MotifOP const & m,
            int parent_index,
            String parent_end_name);

    void
    add_connection(
            int,
            int,
            String const &,
            String const &);


public: //removal functions

    void
    remove_node(int i = -1);

    void
    remove_node_level(int level = -1);


public: //outputing functions

    inline
    void
    to_pdb(
            String const fname = "test.pdb",
            int renumber = -1,
            int close_chains = 0,
            int conect_statements = 0) {

        try {
            _update_merger();
            return merger_->to_pdb(fname, renumber, close_chains, conect_statements);
        }
        catch (MotifMergerException) {
            throw MotifTreeException(
                    "cannot produce merged structure for a pdb it is likely you have created a ring "
                            "with no start, call write_pdbs() to see what the topology would look like");
        }

    }

    String
    topology_to_str();

    String
    to_pretty_str() {
        auto printer = MotifTreePrinter(*this);
        return printer.print_tree(*this);
    }

    String
    to_str();

public: // misc

    void
    _update_merger();

public: //getters

    MotifConnections const &
    connections() { return connections_; }

    util:Beads
    beads();


public: //tree wrappers

    size_t
    size() { return tree_.size(); }

    inline
    data_structure::tree::TreeNodeOP<motif::MotifOP> const &
    get_node(int i) {
        try {
            return tree_.get_node(i);
        }
        catch (data_structure::tree::TreeException) {
            throw MotifTreeException(
                    "cannot get node: " + std::to_string(i) + " in MotifTree it does not exist");
        }
    }

    inline
    data_structure::tree::TreeNodeOP<motif::MotifOP> const &
    get_node(util::Uuid const & uuid) {
        for (auto const & n : tree_) {
            if (n->data()->id() == uuid) {
                return n;
            }
        }
        throw MotifTreeException(
                "cannot get node with uuid no motif has it in this tree");
    }

    inline
    data_structure::tree::TreeNodeOP<motif::MotifOP>
    get_node(String const & m_name) {
        auto node = data_structure::tree::TreeNodeOP<motif::MotifOP>(nullptr);
        for (auto const & n : tree_) {
            if (n->data()->name() == m_name) {
                if (node != nullptr) {
                    throw MotifTreeException(
                            "cannot get node with name: " + m_name + " there is more then one motif "
                                    "with this name");
                }

                node = n;
            }
        }

        if (node == nullptr) {
            throw MotifTreeException(
                    "cannot get node with name: " + m_name + " there is no motif in the tree with "
                            "this name");
        }

        return node;
    }

    inline
    data_structure::tree::TreeNodeOP<motif::MotifOP> const &
    last_node() { return tree_.last_node(); }

    void
    write_pdbs(String const & fname = "nodes");

    inline
    void
    increase_level() { tree_.increase_level(); }

    inline
    void
    decrease_level() { tree_.decrease_level(); }

public: //merger wrappers

    inline
    structure::RNAStructureOP const &
    get_structure() {
        try {
            _update_merger();
            return merger_->get_structure();
        }
        catch (MotifMergerException) {
            throw MotifTreeException(
                    "cannot produce merged structure it is likely you have created a ring with no start"
                            "call write_pdbs() to see what the topology would look like");
        }
    }

    secondary_structure::PoseOP
    secondary_structure() {
        try {
            _update_merger();
            return merger_->secondary_structure();
        }
        catch (MotifMergerException) {
            throw MotifTreeException(
                    "cannot produce merged secondary structure it is likely you have created a ring "
                            "with no start, call write_pdbs() to see what the topology would look like");
        }

    }

    secondary_structure::PoseOP
    designable_secondary_structure() {
        auto ss = merger_->secondary_structure();
        auto ss_r = secondary_structure::ResidueOP(nullptr);

        for (auto const & n : tree_) {
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


public: //option wrappers

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

private: //private option functions
    void
    setup_options();

    void
    update_var_options();

private:
    data_structure::tree::TreeStatic<motif::MotifOP> tree_;
    MotifMergerOP merger_;
    MotifConnections connections_;
    bool sterics_;
    float clash_radius_;
    int update_merger_;
    base::Options options_;
};

typedef std::shared_ptr<MotifTree> MotifTreeOP;

}

#endif /* defined(__RNAMake__motif_tree__) */































