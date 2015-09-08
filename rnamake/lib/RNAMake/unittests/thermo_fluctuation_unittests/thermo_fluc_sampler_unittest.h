//
//  thermo_fluc_sampler_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__thermo_fluc_sampler_unittest__
#define __RNAMake__thermo_fluc_sampler_unittest__

#include <stdio.h>

#include "unittest.h"

class ThermoFlucSamplerUnittest : public Unittest {
public:
    ThermoFlucSamplerUnittest() {}
    
    ~ThermoFlucSamplerUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_sample();

public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__thermo_fluc_sampler_unittest__) */
