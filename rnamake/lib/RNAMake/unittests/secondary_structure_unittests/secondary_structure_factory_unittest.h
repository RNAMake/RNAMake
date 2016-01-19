//
//  secondary_structure_factory_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_factory_unittest__
#define __RNAMake__secondary_structure_factory_unittest__

#include <stdio.h>

#include <stdio.h>

#include "unittest.h"

namespace unittests {


class SecondaryStructureFactoryUnittest : public Unittest {
public:
    SecondaryStructureFactoryUnittest() {}
    
    ~SecondaryStructureFactoryUnittest() {}
    
public:
    
    void
    test_creation();
        
public:
    
    int
    run();
    
    //void
    //run_all();
    
};
    
}



#endif /* defined(__RNAMake__secondary_structure_factory_unittest__) */
