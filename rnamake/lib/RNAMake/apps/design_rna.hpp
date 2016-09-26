//
//  design_rna.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/26/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef design_rna_hpp
#define design_rna_hpp

#include <stdio.h>
#include "base/application.hpp"
#include "motif_state_search/motif_state_search.h"
#include "motif_data_structures/motif_graph.h"

struct EndStateInfo {
    String name;
    int n_pos;
};


class DesignRNAApp : public Application {
public:
    DesignRNAApp() : Application(),
    search_(MotifStateSearch())
    {}
    
    ~DesignRNAApp() {}
    
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
    _setup_from_pdb();
    
    void
    _setup_from_mg();
    
private:
    MotifStateSearch search_;
    MotifGraph mg_;
    EndStateInfo start_, end_;
    
};



#endif /* design_rna_hpp */
