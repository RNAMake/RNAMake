//
//  motif_tree_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_unittest_new__
#define __RNAMake__motif_tree_unittest_new__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace motif_structures {

class MotifTreeUnittest : public Unittest {
public:
    
    MotifTreeUnittest() {}
    
    ~MotifTreeUnittest() {}
    
public:
    
    void
    test_creation();
    
    void
    test_add_motif();
    
    void
    test_remove();
    
public:
    
    int
    run();
    
    //int
    //run_all();
    
private:
    
    
};

}
}


#endif /* defined(__RNAMake__motif_tree_unittest__) */
