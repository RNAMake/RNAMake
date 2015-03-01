//
//  sequence_designer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_designer__
#define __RNAMake__sequence_designer__

#include <stdio.h>
#include "secondary_structure_tree.h"
#include "vienna.h"

class SequenceDesigner {
public:
    SequenceDesigner():
    v_ (Vienna() ) {}
    
    ~SequenceDesigner() {}
    
public:
    
    String
    design(
        String const &,
        String const &);
    
private:
    SecondaryStructureTree ss_tree_;
    Vienna v_;
};

#endif /* defined(__RNAMake__sequence_designer__) */
