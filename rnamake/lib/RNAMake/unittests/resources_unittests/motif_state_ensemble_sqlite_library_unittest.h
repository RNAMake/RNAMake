//
//  motif_state_ensemble_sqlite_library_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_ensemble_sqlite_library_unittest__
#define __RNAMake__motif_state_ensemble_sqlite_library_unittest__

#include <stdio.h>

#include "unittest.h"

class MotifStateEnsembleSqliteLibraryUnittest : public Unittest {
public:
    
    MotifStateEnsembleSqliteLibraryUnittest() {}
    
    ~MotifStateEnsembleSqliteLibraryUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_get();
    
    int
    test_get_random();
    
    int
    test_all();
    
    int
    test_get_multi();
    
    int
    test_contains();
    
public:
    
    int
    run();
    
    
};


#endif /* defined(__RNAMake__motif_state_ensemble_sqlite_library_unittest__) */
