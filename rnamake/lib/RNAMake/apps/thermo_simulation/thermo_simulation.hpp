//
//  thermo_simulation.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef thermo_simulation_hpp
#define thermo_simulation_hpp

#include <stdio.h>

#include <stdio.h>

#include "base/application.hpp"

class ThermoSimulationApp : public Application {
public:
    ThermoSimulationApp() : Application(),
    
    ~ThermoSimulationApp() {}
    
public:
    
    void
    setup_options();
    
    void
    parse_command_line(
                       int,
                       const char **);
    
    void
    run();
    
    
private:
        
};



#endif /* thermo_simulation_hpp */
