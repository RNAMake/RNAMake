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
    SecondaryStructureTree(
        String const & ss,
        String const & seq) {
        String css (ss);
        String cseq (seq);
        _build_tree(css, cseq);

    }
    
    ~SecondaryStructureTree() {
        for(int i = 0; i < nodes_.size(); i++) { delete nodes_[i]; }
    }
    
public:
    SSandSeqOP
    get_ss_and_seq() { return nodes_[0]->get_ss_and_seq(); }
   
    inline
    SecondaryStructureNodes
    get_bulges() const {
        SecondaryStructureNodes bulges;
        for (auto const & n : nodes_) {
            if(n->ss_type() == SSN_TWOWAY) { bulges.push_back(n); }
        }
        return bulges;
    }
    
    
private:
    
    void
    _build_tree(String &, String &);

private:
    SecondaryStructureNodes nodes_;
    
};

#endif /* defined(__RNAMake__secondary_structure_tree__) */
