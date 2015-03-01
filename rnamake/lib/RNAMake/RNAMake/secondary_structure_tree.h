//
//  secondary_structure_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_tree__
#define __RNAMake__secondary_structure_tree__
#include <stdio.h>
#include <iostream>
#include "types.h"
#include "secondary_structure_node.fwd.h"
#include "secondary_structure_node.h"


class SecondaryStructureTree {
public:
    
    SecondaryStructureTree():
    nodes_( SecondaryStructureNodeOPs() )
    {}
    
    SecondaryStructureTree(
        String const & ss,
        String const & seq) {
        String css (ss);
        String cseq (seq);
        _build_tree(css, cseq);

    }
    
    ~SecondaryStructureTree() {}
    
public:
    SSandSeqOP
    get_ss_and_seq() { return nodes_[0]->get_ss_and_seq(); }
   
    inline
    SecondaryStructureNodeOPs
    get_bulges() const {
        SecondaryStructureNodeOPs bulges;
        for (auto const & n : nodes_) {
            if(n->ss_type() == SSN_TWOWAY) { bulges.push_back(n); }
        }
        return bulges;
    }
    
    inline
    SecondaryStructureNodeOPs
    get_designable_bps() const {
        SecondaryStructureNodeOPs bps;
        
        for ( auto const & n : nodes_) {
            if(n->ss_type() == SSN_BP && n->bp_type().compare("NN") == 0) {
                bps.push_back(n);
            }
        }
        
        return bps;
    }
    
    
private:
    
    void
    _build_tree(String &, String &);

private:
    SecondaryStructureNodeOPs nodes_;
    
};

#endif /* defined(__RNAMake__secondary_structure_tree__) */
