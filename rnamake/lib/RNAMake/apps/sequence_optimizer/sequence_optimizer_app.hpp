//
//  sequence_optimizer.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef sequence_optimizer_app_hpp
#define sequence_optimizer_app_hpp

#include <stdio.h>

#include "base/application.hpp"
#include "sequence_optimizer/sequence_optimizer.h"

class SequenceOptimizerApp : public Application {
public:
    SequenceOptimizerApp() : Application(),
    optimizer_(SequenceOptimizer()) {}
    
    ~SequenceOptimizerApp() {}
    
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
    void
    run_from_mg();
    
    
private:
    
    SequenceOptimizer optimizer_;
    
};



#endif /* sequence_optimizer_hpp */
