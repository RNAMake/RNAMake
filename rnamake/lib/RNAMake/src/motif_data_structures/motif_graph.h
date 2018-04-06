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
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_merger.h"

enum MotifGraphStringType {
    OLD,
    MG,
    TOP
};

class MotifGraphException : public std::runtime_error {
public:
    MotifGraphException(
        String const & message) :
    std::runtime_error(message)
    {}
};

class MotifGraph {
private:
    struct _MotifGraphBuildPoint {
        _MotifGraphBuildPoint(
            GraphNodeOP<MotifOP> const & nnode,
            int nend_index ):
        node(nnode),
        end_index(nend_index)
        {}
        
        GraphNodeOP<MotifOP> node;
        int end_index;
    };
    
    typedef std::shared_ptr<_MotifGraphBuildPoint> _MotifGraphBuildPointOP;
    typedef std::vector<_MotifGraphBuildPointOP> _MotifGraphBuildPointOPs;
    
    friend class MotifGraphPrinter;
    
    class MotifGraphPrinter {
    public:
        inline
        MotifGraphPrinter(
            MotifGraph const & mg):
        levels_(std::map<int, int>()),
        node_pos_(std::map<int, int>()),
        nodes_per_level_(std::map<int, int>()),
        branch_length_(25),
        start_pos_(100),
        nodes_(GraphNodeOPs<MotifOP>()){
            
            _setup_node_positions(mg);
        }
        
        ~MotifGraphPrinter() {}
        
    public:
        String
        print_graph(
            MotifGraph const & mg) {
            auto s = String("\n");

            auto nodes_per_level = std::map<int, GraphNodeOPs<MotifOP>>();
            for(auto const & n : mg) {
                nodes_per_level[levels_[n->index()]].push_back(n);
            }
            
            int found = 1, level = 1;
            
            while(found) {
                if(nodes_per_level.find(level) == nodes_per_level.end()) { break; }
                auto node_level = nodes_per_level[level];
                s += _print_level(node_level);
                level += 1;
                
            }
            
            return s;
        }
        
    private:
        String
        _print_level(GraphNodeOPs<MotifOP> const & nodes) {
            auto s = String("");
            auto strings = std::vector<Strings>();
            
            //sort nodes by pos, using a lambda statement and std::pair
            //probably can be cleaned up
            auto pairs = std::vector<std::pair<GraphNodeOP<MotifOP>, int>>();
            for(auto const & n : nodes) {
                auto p = std::pair<GraphNodeOP<MotifOP>, int>(n, node_pos_[n->index()]);
                pairs.push_back(p);
            }
            
            std::sort(pairs.begin(), pairs.end(),
                      [](std::pair<GraphNodeOP<MotifOP>, int> const & p1,
                         std::pair<GraphNodeOP<MotifOP>, int> const & p2) {
                          return p1.second < p2.second; });
            
            
            auto sorted_nodes = GraphNodeOPs<MotifOP>();
            for(auto const & p : pairs) { sorted_nodes.push_back(p.first); }

            
            for(auto const & n : sorted_nodes) {
                auto strs = Strings();
                if(n->parent() != nullptr) {
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
            for(int i = 0; i < strings[0].size(); i++) {
                auto transposed = Strings();
                for(auto const & strs: strings) {
                    transposed.push_back(strs[i]);
                }
                transposed_strings.push_back(transposed);
            }
            
            for(auto const & strs : transposed_strings) {
                int current_pos = 0;
                int j = 0;
                for(auto const & n : sorted_nodes) {
                    int pos = node_pos_[n->index()];
                    int diff = pos - current_pos;
                    for(int i = 0; i < diff; i++) { s += " "; }
                    s += strs[j];
                    current_pos = pos + strs[j].length();
                    j += 1;
                }
                s += "\n";
            }
            
            int current_pos = 0;
            int hit = 0;
            for(auto const & n : sorted_nodes) {
                auto children = GraphNodeOPs<MotifOP>();
                int j = -1;
                for(auto const & c : n->connections()) {
                    j++;
                    if(j == 0) { continue; }
                    if(c != nullptr) {
                        children.push_back(c->partner(n->index()));
                    }
                }
                
                if(children.size() > 1) {
                    hit = 1;
                    auto min = node_pos_[children[0]->index()];
                    auto max = node_pos_[children.back()->index()];
                    auto diff = min+1 - current_pos;
                    for(int i = 0; i < diff; i++) { s += " "; }
                    for(int i = 0; i < max-min-1; i++) { s += "_"; }
                    current_pos = max;
                }
            }
            
            if(hit) { s += "\n"; }
            
            return s;
        }
        
        void
        _setup_node_positions(
            MotifGraph const & mg) {
            
            _assign_node_levels(mg);
            
            for(auto const & n : nodes_) {
                auto n_level = levels_[n->index()];
                if(nodes_per_level_.find(n_level) == nodes_per_level_.end()) {
                    nodes_per_level_[n_level] = 0;
                }
                nodes_per_level_[n_level] += 1;
            }
            
            int i = -1;
            for(auto const & n : nodes_) {
                i++;
                if(i == 0) {
                    node_pos_[n->index()] = start_pos_;
                }
                
                auto children = GraphNodeOPs<MotifOP>();
                int j = -1;
                for(auto const & c : n->connections()) {
                    j++;
                    if(j == 0) { continue; }
                    if(c != nullptr) {
                        auto child = c->partner(n->index());
                        if(child->connections()[0] != c) { continue; }
                        if(node_pos_.find(child->index()) != node_pos_.end()) { continue; }
                        
                        children.push_back(child);
                    }
                }
                
                if(children.size() == 1) {
                    node_pos_[children[0]->index()] = node_pos_[n->index()];
                }
                else if(children.size() == 2) {
                    auto level = levels_[n->index()];
                    auto nodes_per_level = nodes_per_level_[level+1];
                    auto extra = nodes_per_level - 2;
                    auto parent_pos = node_pos_[n->index()];
                    if(extra == 0) {
                        node_pos_[children[0]->index()] = parent_pos - branch_length_;
                        node_pos_[children[1]->index()] = parent_pos + branch_length_;
                    }
                    else {
                        node_pos_[children[0]->index()] = parent_pos - branch_length_ / extra;
                        node_pos_[children[1]->index()] = parent_pos + branch_length_ / extra;
                    }
                }
                else if(children.size() > 2) {
                    throw MotifGraphException(
                        "Greater then two children is not supported for pretty_printing");
                }
                
            }
        }
       
        void
        _assign_node_levels(
            MotifGraph const & mg) {
            
            auto nodes = mg.unaligned_nodes();
            
            if(nodes.size() == 0) {
                throw MotifGraphException(
                    "cannot find a place to start printing in motif_graph to_pretty_str");
            }
            
            auto start = nodes[0]->index();
            
            auto n = GraphNodeOP<MotifOP>();
            for(auto it = mg.graph_.transverse(mg.graph_.get_node(start));
                it != mg.graph_.end();
                ++it) {
                
                n = (*it);
                
                if(levels_.size() == 0) {
                    levels_[n->index()] = 1;
                }
                else {
                    auto c = n->connections()[0];
                    auto parent = c->partner(n->index());
                    levels_[n->index()] = levels_[parent->index()] + 1;
                }
                nodes_.push_back(n);

            }
        }
        
        

    private:
        std::map<int, int> levels_;
        std::map<int, int> node_pos_;
        std::map<int, int> nodes_per_level_;
        int branch_length_;
        int start_pos_;
        GraphNodeOPs<MotifOP> nodes_;
        
        
    };

    
public: //construction
    
    MotifGraph();
    
    MotifGraph(
        String const &);
    
    MotifGraph(
        String const &,
        MotifGraphStringType const &);
    
    MotifGraph(
        MotifGraph const &);

    ~MotifGraph() { }

public: //setup helpers
    
    void
    update_indexes(std::map<int, int> const &);
    
public: //iterators
    
    typedef typename GraphStatic<MotifOP>::iterator iterator;
    typedef typename GraphStatic<MotifOP>::const_iterator const_iterator;
    
    iterator begin() { return graph_.begin(); }
    iterator end()   { return graph_.end(); }
    
    const_iterator begin() const { return graph_.begin(); }
    const_iterator end()   const { return graph_.end(); }
    
private://add function helpers
    
    GraphNodeOP<MotifOP>
    _get_parent(
        String const &,
        int);
    
    int
    _add_motif_to_graph(
        MotifOP &,
        GraphNodeOP<MotifOP> const &,
        int);

    Ints
    _get_available_parent_end_pos(
        GraphNodeOP<MotifOP> const &,
        int);
    
    void
    _add_motif_tree(
        MotifTreeOP const &,
        int,
        int);
    
    int
    _steric_clash(
        MotifOP const &);
    
    
    int
    _get_connection_end(
        GraphNodeOP<MotifOP> const &,
        String const &);
    
    void
    _align_motifs_all_motifs();

public: //add functions
    
    int
    add_motif(
        MotifOP const & m,
        int parent_index = -1,
        int parent_end_index = -1,
        int orphan = 0);
    
    int
    add_motif(
        MotifOP const &,
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
    GraphNodeOP<MotifOP> const &
    last_node() { return graph_.last_node(); }
    
    inline
    GraphNodeOP<MotifOP>
    oldest_node() { return graph_.oldest_node(); }
    
    inline
    GraphNodeOP<MotifOP> const &
    get_node(int i) const { return graph_.get_node(i); }
    
    inline
    GraphNodeOP<MotifOP> const
    get_node(Uuid const & uuid) const {
        for(auto const & n : graph_) {
            if(n->data()->id() == uuid) {
                return n;
            }
        }
        throw MotifGraphException("cannot get node with uuid no motif has it in this graph");
    }
    
    inline
    GraphNodeOP<MotifOP> const
    get_node(String const & m_name) const {
        auto node = GraphNodeOP<MotifOP>(nullptr);
        for(auto const & n : graph_) {
            if(n->data()->name() == m_name) {
                if(node != nullptr) {
                    throw MotifGraphException(
                        "cannot get node with name: " + m_name + " there is more then one motif "
                        "with this name");
                }
                
                node = n;
            }
        }
        
        if(node == nullptr) {
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
    replace_helical_sequence(sstruct::PoseOP const &);
    
    inline
    void
    replace_helical_sequence(String const & seq) {
        auto dss = designable_secondary_structure();
        dss->replace_sequence(seq);
        replace_helical_sequence(dss);
    }
    
    sstruct::PoseOP
    designable_secondary_structure() {
        _update_merger();
        auto ss = merger_->secondary_structure();
        auto ss_r = sstruct::ResidueOP(nullptr);
        
        for(auto const & n : graph_) {
            if(n->data()->name() != "HELIX.IDEAL") { continue; }
            
            for(auto const & r : n->data()->residues()) {
                ss_r= ss->get_residue(r->uuid());
                if(ss_r != nullptr) {
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
    
    String
    to_pretty_str() {
        auto printer = MotifGraphPrinter(*this);
        return printer.print_graph(*this);
    }
    
public: // misc functions
    
    void
    _update_align_list();
    
    void
    _update_merger();

public: // getters
    
    inline
    Beads
    beads() {
        Beads beads;
        for(auto const & n : graph_.nodes()) {
            std::copy(n->data()->beads().begin(),
                      n->data()->beads().end(),
                      std::inserter(beads, beads.end()));
        }
        return beads;
    }
    
    BasepairOP const &
    get_available_end(int);
    
    BasepairOP const &
    get_available_end(
        int,
        String const &);
    
    BasepairOP 
    get_available_end(
        String const &,
        String const &);
    
    GraphNodeOPs<MotifOP> const
    unaligned_nodes() const;
    
    _MotifGraphBuildPointOPs
    get_build_points();
    
    // this is bad
    std::map<int, int>
    aligned() { return aligned_; }
    
    GraphConnectionOPs<MotifOP> const &
    connections() { return graph_.connections(); }
    
    
    
public: //Motif Merger Wrappers
    
    inline
    RNAStructureOP const &
    get_structure() {
        try {
            _update_merger();
            return merger_->get_structure();
        }
        catch(MotifMergerException) {
            throw MotifGraphException(
                "cannot produce merged structure it is likely you have created a ring with no start"
                "call write_pdbs() to see what the topology would look like");
        }
    }
    
    sstruct::PoseOP
    secondary_structure() {
        try {
            _update_merger();
            return merger_->secondary_structure();
        }
        catch(MotifMergerException) {
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
        catch(MotifMergerException) {
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
    GraphStatic<MotifOP> graph_;
    MotifMergerOP merger_;
    Options options_;
    std::map<int, int> aligned_;
    GraphNodeOPs<MotifOP> align_list_;
    int update_merger_;
    int update_align_list_;
    //options
    float clash_radius_;
    bool sterics_;

};

typedef std::shared_ptr<MotifGraph> MotifGraphOP;

#endif /* defined(__RNAMake__motif_graph__) */
