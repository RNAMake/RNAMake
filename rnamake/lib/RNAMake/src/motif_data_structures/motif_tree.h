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
#include "motif_data_structures/motif_merger.h"
#include "motif_data_structures/motif_connection.h"

class MotifTreeException : public std::runtime_error {
public:
    MotifTreeException(
        String const & message) :
    std::runtime_error(message)
    {}
};


class MotifTree  {
private:
    class MotifTreePrinter {
    public:
        inline
        MotifTreePrinter():
        levels_(std::map<int, int>()),
        node_pos_(std::map<int, int>())
        {}
        
        ~MotifTreePrinter() {}
        
    public:
        String
        print_tree(
            MotifTree const & mt) {
            auto s = String("");
            
            _assign_node_levels(mt);
            return s;
            
            int c_level = 1;
            while(1) {
                auto node_level = TreeNodeOPs<MotifOP>();
                for(auto const & n : mt) {
                    if(levels_[n->index()] == c_level) {
                        node_level.push_back(n);
                    }
                }
                
                if(c_level == 1) {                }
                
                if(node_level.size() == 0) { break; }
                
                //_print_node_level(node_level);
                
                c_level++;
    
                
            }
            
            return s;
        }
        
        
    private:
        void
        _assign_node_levels(
            MotifTree const & mt) {
            
            max_level_ = 0;
            
            for(auto const & n : mt) {
                if(levels_.size() == 0) {
                    levels_[n->index()] = 1;
                    continue;
                }
                else {
                    auto parent_level = levels_[n->parent_index()];
                    levels_[n->index()] = parent_level + 1;
                    if(parent_level + 1 > max_level_) {
                        max_level_ = parent_level + 1;
                    }
                }
            }
            
            int c_level = 1;
            int max_node_level = 0;
            while(1) {
                int node_level = 0;
                for(auto const & n : mt) {
                    if(levels_[n->index()] == c_level) { node_level++; }
                }
                
                if(max_node_level < node_level) { max_node_level = node_level; }
                
                c_level++;

                if(node_level == 0) { break; }
                
            }
            
            
            int branchLen = 10*max_node_level;
            int startLen = 20*max_node_level;
            c_level = 1;
            
            int found = 1;
            
            while(found) {
                found = 0;
                for(auto const & n : mt) {
                    if(levels_[n->index()] == c_level) {
                        found = 1;
                        if(c_level == 1) {
                            node_pos_[n->index()] = startLen;
                        }
                        auto children = TreeNodeOPs<MotifOP>();
                        for(auto const & c : n->children()) {
                            if(c != nullptr) { children.push_back(c); }
                        }

                        if(children.size() == 1) {
                            
                        }
                    
                    }
                }
                
                c_level++;
                
                
            }
            
            
        }
        
        void
        _print_node_level(
            TreeNodeOPs<MotifOP> const & nodes) {
            
            auto all_children = std::vector<TreeNodeOPs<MotifOP>>();
           
            for(auto const & n : nodes) {
                auto children = TreeNodeOPs<MotifOP>();
                for(auto const & c : n->children()) {
                    if(c != nullptr) { children.push_back(c); }
                }
                all_children.push_back(children);
            }
            
            //have to do one line at a time for all nodes
            //print index and name of node
            for(int i = 0; i < nodes.size(); i++ ) {
                auto pos = node_pos_[nodes[i]->index()];
                if(all_children[i].size() < 2) {
                    std::cout << std::setw(pos);
                    std::cout << nodes[i]->index() << " - " << nodes[i]->data()->name();
                }
                else {
                    std::cout << std::setw(pos-21) << "";
                    std::cout << std::setfill('_') << std::setw(20) << "";
                    std::cout << nodes[i]->index() << " - " << nodes[i]->data()->name();
                    std::cout << std::setfill('_') << std::setw(20) << "";
                    std::cout << std::setfill(' ');
                }
                
            }
            std::cout << std::endl;
     
            //print end_name of 0 index
            for(int i = 0; i < nodes.size(); i++ ) {
                auto pos = node_pos_[nodes[i]->index()];
                if(all_children[i].size() < 2) {
                    std::cout << std::setw(pos) << "";
                    std::cout <<" - " << nodes[i]->data()->ends()[0]->name();
                }
                if(all_children[i].size() == 2) {
                    std::cout << std::setw(pos-20) << "|";
                    std::cout << std::setw(pos+20) << "|";
                }
                
            }
            
            std::cout << std::endl;
            
            //print first connection line
            for(int i = 0; i < nodes.size(); i++ ) {
                auto pos = node_pos_[nodes[i]->index()];
                if(all_children[i].size() == 1) {
                    std::cout << std::setw(pos) <<"|";
                }
                if(all_children[i].size() == 2) {
                    std::cout << std::setw(pos-20) << "|";
                    std::cout << std::setw(pos+20) << "|";
                }
                
            }
            
            std::cout << std::endl;
            
            //print first second line
            for(int i = 0; i < nodes.size(); i++ ) {
                auto pos = node_pos_[nodes[i]->index()];
                if(all_children[i].size() == 1) {
                    std::cout << std::setw(pos-1) << "";
                    std::cout << "| - ";
                    std::cout << all_children[i][0]->parent_end_index() << " - ";
                    std::cout << nodes[i]->data()->end_name(all_children[i][0]->parent_end_index());
                    
                }
                
            }
            
            std::cout << std::endl;

            for(auto const & n : nodes) {
                auto children = TreeNodeOPs<MotifOP>();
                for(auto const & c : n->children()) {
                    if(c != nullptr) {
                        children.push_back(c);
                    }
                }
                
                auto pos = node_pos_[n->index()];
                if(children.size() == 1) {
                    std::cout << std::setw(pos);
                    std::cout << "|";
                    node_pos_[children[0]->index()] = pos;
                }
            }
            std::cout << std::endl;

        }
        
    private:
        int max_level_;
        std::map<int, int> levels_;
        std::map<int, int> node_pos_;
    };
    
public:
    
    MotifTree():
    tree_(TreeStatic<MotifOP>()),
    merger_(MotifMerger()),
    connections_(MotifConnections()),
    options_(Options())    
    { setup_options(); }
    
    MotifTree(
        String const &);
    
    ~MotifTree() {}
    
public: //iterators
    
    typedef typename TreeStatic<MotifOP>::iterator iterator;
    typedef typename TreeStatic<MotifOP>::const_iterator const_iterator;
    
    iterator begin() { return tree_.begin(); }
    iterator end()   { return tree_.end(); }
    
    const_iterator begin() const { return tree_.begin(); }
    const_iterator end()   const { return tree_.end(); }
    
public:
    String
    to_pretty_str();
    
    
public: //tree wrappers
    
    size_t
    size() { return tree_.size(); }
    
    inline
    TreeNodeOP<MotifOP> const &
    get_node(int i) {
        try {
            return tree_.get_node(i);
        }
        catch(TreeException) {
            throw MotifTreeException(
                "cannot get node: " + std::to_string(i) + " in MotifTree it does not exist");
        }
    }
    
    inline
    TreeNodeOP<MotifOP> const &
    get_node(Uuid const & uuid) {
        for(auto const & n : tree_) {
            if(n->data()->id() == uuid) {
                return n;
            }
        }
        throw MotifTreeException(
            "cannot get node with uuid no motif has it in this tree");
    }
    
    inline
    TreeNodeOP<MotifOP> 
    get_node(String const & m_name) {
        auto node = TreeNodeOP<MotifOP>(nullptr);
        for(auto const & n : tree_) {
            if(n->data()->name() == m_name) {
                if(node != nullptr) {
                    throw MotifTreeException(
                        "cannot get node with name: " + m_name + " there is more then one motif "
                        "with this name");
                }
                
                node = n;
            }
        }
        
        if(node == nullptr) {
            throw MotifTreeException(
                "cannot get node with name: " + m_name + " there is no motif in the tree with "
                "this name");
        }
        
        return node;
    }
    
    inline
    TreeNodeOP<MotifOP> const &
    last_node() { return tree_.last_node(); }
    
    void
    write_pdbs(String const & fname = "nodes");
    
public: //merger wrappers
    
    inline
    RNAStructureOP const &
    get_structure() {
        try { return merger_.get_structure(); }
        catch(MotifMergerException) {
            throw MotifTreeException(
                "cannot produce merged structure it is likely you have created a ring with no start"
                "call write_pdbs() to see what the topology would look like");
        }
    }
    
    sstruct::PoseOP
    secondary_structure() {
        try { return merger_.secondary_structure(); }
        catch(MotifMergerException) {
            throw MotifTreeException(
                "cannot produce merged secondary structure it is likely you have created a ring "
                "with no start, call write_pdbs() to see what the topology would look like");
        }

    }
    
    inline
    void
    to_pdb(
        String const fname = "test.pdb",
        int renumber = -1) {
        
        try { return merger_.to_pdb(fname, renumber); }
        catch(MotifMergerException) {
            throw MotifTreeException(
                "cannot produce merged structure for a pdb it is likely you have created a ring "
                "with no start, call write_pdbs() to see what the topology would look like");
        }
        
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

    
public: //add motif interface
    
    int
    add_motif(
        MotifOP const & m,
        int parent_index = -1,
        int parent_end_index = -1);
    
    int
    add_motif(
        MotifOP const & m,
        int parent_index,
        String parent_end_name);
    

private: // add motif helper functions
    
    TreeNodeOP<MotifOP>
    _get_parent(
        String const & m_name,
        int parent_index) {
        
        auto parent = tree_.last_node();
        
        //catch non existant parent
        try {
            if(parent_index != -1) { parent = tree_.get_node(parent_index); }
        }
        catch(TreeException e) {
            throw MotifTreeException("could not add motif: " + m_name + " with parent: "
                                     + std::to_string(parent_index) + "there is no node with" +
                                     "that index");
        }
        
        return parent;
    }
    
    Ints
    _get_available_parent_end_pos(
        String const & m_name,
        TreeNodeOP<MotifOP> const & parent,
        int parent_end_index) {
        
        auto avail_pos = Ints();
        
        if(parent_end_index != -1) {
            int avail = parent->available_pos(parent_end_index);
            if(!avail) {
                throw MotifTreeException(
                    "could not add motif: " + m_name+ " with parent: " + std::to_string(parent->index()) +
                    " since the parent_end_index supplied " + std::to_string(parent_end_index) +
                                         " is filled");
            }
            else {
                auto name = parent->data()->ends()[parent_end_index]->name();
                if(connections_.in_connection(parent->index(), name)) {
                    throw MotifTreeException(
                        "could not add motif: " + m_name + " with parent: " + std::to_string(parent->index()) +
                        " since the parent_end_index supplied " + std::to_string(parent_end_index) +
                        " is in a connection");
                }
                
                avail_pos.push_back(parent_end_index);
            }
        }
        
        else {
            avail_pos = parent->available_children_pos();
            if(avail_pos.size() == 1 && avail_pos[0] == parent->data()->block_end_add()) {
                throw MotifTreeException(
                    "could not add motif: " + m_name + " with parent: " + std::to_string(parent->index()) +
                    " since it has no free ends to add too");
            }
            
        }
        
        return avail_pos;

    }
    
    int
    _steric_clash(
        MotifOP const &);

public:
    
    void
    add_connection(
        int,
        int,
        String const &,
        String const &);
    
    String
    topology_to_str();

    void
    remove_node(
        int i=-1) {
        
        if(i == -1) {
            i = last_node()->index();
        }
        
        try {
            auto n = get_node(i);
            tree_.remove_node(n);
            merger_.remove_motif(n->data());
            connections_.remove_connections_to(i);
        
        }
        catch(MotifTreeException) {
            throw MotifTreeException(
                "cannot remove node with index: " + std::to_string(i) + " as it does not exist");
        }
    }
    
private:
    void
    setup_options();
    
    void
    update_var_options();
    
private:
    TreeStatic<MotifOP> tree_;
    MotifMerger merger_;
    MotifConnections connections_;
    bool sterics_;
    float clash_radius_;
    Options options_;
};

typedef std::shared_ptr<MotifTree> MotifTreeOP;


#endif /* defined(__RNAMake__motif_tree__) */































