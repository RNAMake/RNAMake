//
//  motif_graph_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_graph_unittest__
#define __RNAMake__motif_graph_unittest__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace motif_structures {

class MotifGraphUnittest : public Unittest {
public:
    
    MotifGraphUnittest() {}
    
    ~MotifGraphUnittest() {}
    
public:
    
    void
    test_creation();
    
    void
    test_add_motif();
    
    void
    test_remove();
    
    void
    test_remove_2();
    
    void
    test_copy();
    
    void
    test_replace_ideal_helices();
    
    void
    test_secondary_structure();
    
public:
    
    int
    run();
    
    //int
    //run_all();
    
private:
    
    
};

}
}



#endif /* defined(__RNAMake__motif_graph_unittest__) */
