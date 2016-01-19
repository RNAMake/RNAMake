//
//  motif_topology_unittests.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_topology_unittests__
#define __RNAMake__motif_topology_unittests__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace motif_structures {

class MotifTopologyUnittest : public Unittest {
public:
    
    MotifTopologyUnittest() {}
    
    ~MotifTopologyUnittest() {}
    
    int
    size() { return 9; }
    
public:
    
    void
    test_to_tree();
    
    void
    test_to_tree_complex();
    
public:
    
    int
    run();
    
    int
    run_all();
    
private:
    
    
};

}
}



#endif /* defined(__RNAMake__motif_topology_unittests__) */
