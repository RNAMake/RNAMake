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
#include "secondary_structure/secondary_structure_parser.h"
#include "secondary_structure/pose.h"

namespace sstruct {

class SecondaryStructureFactory {
public:
    SecondaryStructureFactory():
    parser_(SecondaryStructureParser())
    {}
    
    ~SecondaryStructureFactory() {}
    
public:
    
    inline
    MotifOP
    motif(
        String const & sequence,
        String const & dot_bracket) {
        return parser_.parse_to_motif(sequence, dot_bracket);
    }
    
    inline
    PoseOP
    pose(
        String const & sequence,
        String const & dot_bracket) {
        return parser_.parse_to_pose(sequence, dot_bracket);
    }
    
    /*SecondaryStructureOP
    get_structure(
        String const & sequence,
        String const & dot_bracket) {
        
        auto sstree = SS_Tree(sequence, dot_bracket);
        auto ss = std::make_shared<SecondaryStructure>(sstree.secondary_structure());
        _get_basepairs(sstree, ss);
        _get_motifs(sstree, ss);
        
        return ss;
        
    }*/

private:
    SecondaryStructureParser parser_;
    


};

}

#endif /* defined(__RNAMake__secondary_structure_factory__) */
