//
//  thermo_fluc_simulation_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__thermo_fluc_simulation_unittest__
#define __RNAMake__thermo_fluc_simulation_unittest__

#include <stdio.h>


#include "unittest.h"

class ThermoFlucSimulationUnittest : public Unittest {
public:
    ThermoFlucSimulationUnittest() {}
    
    ~ThermoFlucSimulationUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_run();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__thermo_fluc_simulation_unittest__) */
