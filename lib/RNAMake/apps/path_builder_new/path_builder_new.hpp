//
//  path_builder_new.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef path_builder_new_hpp
#define path_builder_new_hpp

#include <stdio.h>

#include "base/application.hpp"
#include "motif_search/motif_state_search.h"
#include "motif_data_structure/motif_graph.h"

struct EndStateInfo {
    String name;
    int n_pos;
};


class PathBuilderNewApp : public base::Application {
public:
    PathBuilderNewApp() : base::Application(),
    search_(motif_search::MotifStateSearch())
    {}
    
    ~PathBuilderNewApp() {}
    
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
    _setup_from_motif();
    
    void
    _setup_from_mg();
    
private:
    motif_search::MotifStateSearch search_;
    motif_data_structure::MotifGraph mg_;
    EndStateInfo start_, end_;
    
};




#endif /* path_builder_new_hpp */
