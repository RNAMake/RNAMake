//
//  segmenter_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__segmenter_unittest__
#define __RNAMake__segmenter_unittest__

#include <stdio.h>

#include "unittest.h"

class SegmenterUnittest : public Unittest {
public:
    
    SegmenterUnittest() {}
    
    ~SegmenterUnittest() {}
    
public:
    
    int
    test_creation();
    
public:
    
    int
    run();
    
    
};



#endif /* defined(__RNAMake__segmenter_unittest__) */
