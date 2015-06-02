//
//  motif_tree_topology.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_topology__
#define __RNAMake__motif_tree_topology__

#include <stdio.h>
#include "base/types.h"
#include "secondary_structure/ss_tree.h"
#include "data_structure/tree/tree.h"


struct MotifTopology {
public:
    enum class MotifType { BASEPAIR, TWOWAY, NWAY, HAIRPIN, TCONTACT};
    MotifTopology(
        MotifType const & ntype,
        String const & nseq,
        String const & nss):
    type(ntype),
    seq(nseq),
    ss(nss)
    {}

    MotifType type;
    String seq;
    String ss;
};

class MotifTreeTopology {
public:
    MotifTreeTopology(SS_Tree const &);
    
private:
    String
    combine_seq(SS_Nodes const &);
   
private:
    Tree<MotifTopology*> tree_;
    
};
 

#endif /* defined(__RNAMake__motif_tree_topology__) */
