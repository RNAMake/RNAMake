//
//  secondary_structure_factory.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_factory__
#define __RNAMake__secondary_structure_factory__

#include <stdio.h>

//RNAMake Headers
#include "secondary_structure/secondary_structure.h"
#include "secondary_structure/ss_tree.h"

namespace sstruct {

class SecondaryStructureFactory {
public:
    SecondaryStructureFactory() {}
    
    ~SecondaryStructureFactory() {}
    
public:
    
    SecondaryStructureOP
    get_structure(
        String const & sequence,
        String const & dot_bracket) {
        
        auto sstree = SS_Tree(sequence, dot_bracket);
        auto ss = std::make_shared<SecondaryStructure>(sstree.secondary_structure());
        _get_basepairs(sstree, ss);
        _get_motifs(sstree, ss);
        
        return ss;
        
    }

private:
    
    void
    _get_basepairs(
        SS_Tree const &,
        SecondaryStructureOP &);
    
    void
    _get_motifs(
        SS_Tree const &,
        SecondaryStructureOP &);
    
    ChainOPs
    _get_chains(
        SecondaryStructureOP const &,
        std::vector<TreeNodeOP<SS_NodeDataOP>> const &);

};

}

#endif /* defined(__RNAMake__secondary_structure_factory__) */
