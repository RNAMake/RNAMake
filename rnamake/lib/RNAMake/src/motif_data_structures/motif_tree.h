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
#include "base/base.h"
#include "data_structure/tree/tree.h"
#include "motif/motif.h"
#include "motif_data_structures/motif_merger.h"

class MotifTreeException : public std::runtime_error {
public:
    MotifTreeException(
        String const & message) :
    std::runtime_error(message)
    {}
};


class MotifTree : public Base  {
public:
    
    MotifTree():
    tree_(TreeStatic<MotifOP>()),
    merger_(MotifMerger()),
    clash_radius_(2.5),
    sterics_(1)
    {}
    
public: //iterators
    
    typedef typename TreeStatic<MotifOP>::iterator iterator;
    typedef typename TreeStatic<MotifOP>::const_iterator const_iterator;
    
    iterator begin() { return tree_.begin(); }
    iterator end()   { return tree_.end(); }
    
    const_iterator begin() const { return tree_.begin(); }
    const_iterator end()   const { return tree_.end(); }
    
public:
    
    size_t
    size() { return tree_.size(); }
    
    inline
    RNAStructureOP const &
    get_structure() {
        return merger_.get_structure();
    }
    
    inline
    TreeNodeOP<MotifOP> const &
    get_node(int i) { return tree_.get_node(i); }
    
    void
    write_pdbs(String const & fname = "nodes");

public: //add motif interface
    
    int
    add_motif(
        MotifOP const & m,
        int parent_index = -1,
        int parent_end_index = -1);
    
    int
    add_motif(
        String const & m_name,
        int parent_index = -1,
        int parent_end_index = 1);
    
    int
    add_motif(
        String const & m_name,
        String const & m_end_name,
        int parent_index = -1,
        int parent_end_index = 1);

public:
    
    int
    _steric_clash(
        MotifOP const &);
    
private:
    void
    setup_options();
    
    void
    update_var_options();
    
private:
    TreeStatic<MotifOP> tree_;
    MotifMerger merger_;
    int clash_radius_, sterics_; 
};

typedef std::shared_ptr<MotifTree> MotifTreeOP;


#endif /* defined(__RNAMake__motif_tree__) */































