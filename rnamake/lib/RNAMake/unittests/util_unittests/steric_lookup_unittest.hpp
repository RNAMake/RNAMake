//
//  steric_lookup_unittest.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef steric_lookup_unittest_hpp
#define steric_lookup_unittest_hpp

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"
#include "util/random_number_generator.h"

class StericLookupUnittest : public Unittest {
public:
    StericLookupUnittest():
    rng_(RandomNumberGenerator())
    {}
    
    
    ~StericLookupUnittest() {}
    
    
public:
    
    int
    test_creation();
    
    int
    test_add_point();
    
    int
    test_add_point_2();
    
    int
    test_add_point_3();

public:
    
    int
    run();
    
    int
    run_all();
    
private:
    
    Point const &
    rand_point(int scale = 10);
    
    
private:
    RandomNumberGenerator rng_;
    Point test_p_;

    
};


#endif /* steric_lookup_unittest_hpp */
