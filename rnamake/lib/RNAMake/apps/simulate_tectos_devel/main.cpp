//
//  main.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 10/31/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>


#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "simulate_tectos_devel.h"


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);
    
    //load extra motifs being used
    String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path+"GAAA_tetraloop");
    resources::Manager::instance().add_motif(base_path+"GGAA_tetraloop");
    
    
    auto app = SimulateTectosApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;
}

