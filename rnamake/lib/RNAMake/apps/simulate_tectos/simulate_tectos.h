//
//  simulate_tectos.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__simulate_tectos__
#define __RNAMake__simulate_tectos__

#include <stdio.h>

//RNAMake Headers
#include "base/base.h"

class SimulateTectos {
public:
    SimulateTectos(
        String const &,
        String const &,
        String const &,
        String const &);
    
public:
    
    void
    run();
    
};

#endif /* defined(__RNAMake__simulate_tectos__) */
