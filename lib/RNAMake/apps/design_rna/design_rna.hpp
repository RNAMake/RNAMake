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
#include "util/steric_lookup.hpp"
#include "motif_search/motif_state_search.h"
#include "motif_data_structure/motif_graph.h"

class DesignRNAAppException : public std::runtime_error {
public:
    DesignRNAAppException(
        String const & message):
    std::runtime_error(message)
    {}
};


struct EndStateInfo {
    String name;
    int n_pos;
};


class DesignRNAApp : public base::Application {
public:

    DesignRNAApp();
    
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
    _setup_sterics();
    
    void
    _setup_from_pdb();
    
    void
    _setup_from_mg();

    std::shared_ptr<MSS_Path>
    _setup_path();

    std::vector<motif::MotifStateOPs>
    _get_libraries();

private:
    motif_search::MotifStateSearch search_;
    motif_data_structure::MotifGraphOP mg_;
    EndStateInfo start_, end_;
    util::StericLookup lookup_;
    sequence_optimization::SequenceOptimizer3D optimizer_;
    
    
};



#endif /* design_rna_hpp */
