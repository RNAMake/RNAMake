//
//  sequence_optimization.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef sequence_optimizer_app_hpp
#define sequence_optimizer_app_hpp

#include <stdio.h>

#include "base/application.hpp"
#include "sequence_optimization/sequence_optimizer_3d.hpp"

class SequenceOptimizerAppException : public std::runtime_error {
public:
    SequenceOptimizerAppException(
            String const & message):
            std::runtime_error(message)
    {}
};


struct NodeIndexandEdge {
    int ni; //node index
    int ei; //end index
};

struct ConnectionTemplate {
    NodeIndexandEdge start;
    NodeIndexandEdge end;
    String type;
};

class SequenceOptimizerApp : public base::Application {
public:
    SequenceOptimizerApp() : base::Application(),
    optimizer_(sequence_optimization::SequenceOptimizer3D()) {}
    
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
    _get_end_connections(
            motif_data_structure::MotifGraphOP);

    ConnectionTemplate
    _parse_end_commandline_args();

    sequence_optimization::SequenceOptimizerScorerOP
    _setup_optimizer_scorer();


    
    
private:
    std::vector<ConnectionTemplate> connections_;
    sequence_optimization::SequenceOptimizer3D optimizer_;
    
    
};



#endif /* sequence_optimizer_hpp */
