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
        int parent_index = -1) {
        
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
    
    inline
    TreeNodeOP<MotifOP>  const &
    get_node_by_id(
        Uuid const & uuid) {
        for(auto const & n : tree_) {
            if(n->data()->id() == uuid) {
                return n;
            }
        }
        throw MotifTreeException("could not find node by id");
    }
    
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































