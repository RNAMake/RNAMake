//
//  motif_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_unittest__
#define __RNAMake__motif_unittest__

#include <stdio.h>

#include "unittest.h"
#include "motif/motif.h"

class MotifUnittest : public Unittest {
public:
    
    MotifUnittest();
    
    ~MotifUnittest() {}
    
public:
    
    int
    test_copy();
    
    int
    test_to_str();
    
    int
    test_secondary_structure();
    
    int
    test_creation_from_dir();
    
    int
    test_get_basepair_by_name();
    
public:
    
    int
    run();
    
private:
    
    Motif m_;
    
    
};



#endif /* defined(__RNAMake__motif_unittest__) */
