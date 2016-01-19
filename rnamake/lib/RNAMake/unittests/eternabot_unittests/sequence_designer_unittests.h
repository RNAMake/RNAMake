//
//  sequence_designer_unittests.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_designer_unittests__
#define __RNAMake__sequence_designer_unittests__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace eternabot_unittests {

class SequenceDesignerUnittest : public Unittest {
public:
    
    SequenceDesignerUnittest() {}
    
    ~SequenceDesignerUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_design();

public:
    
    int
    run();
    
    int
    run_all();
    
};


    }
}

#endif /* defined(__RNAMake__sequence_designer_unittests__) */
