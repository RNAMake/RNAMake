//
//  mini_ttr.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/16/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__mini_ttr__
#define __RNAMake__mini_ttr__

#include <stdio.h>

//RNAMake Headers
#include "base/option.h"
#include "motif_state_search/motif_state_search.h"


Options
parse_command_line(
    int argc,
    const char ** argv);

class MiniTTRPathFollow {
public:
    
    MiniTTRPathFollow():
    search_(MotifStateSearch()) {
        
    }

    ~MiniTTRPathFollow() {}
    
public:
    void
    setup(Options & opts) {
        search_.set_option_value("max_node_level", 40);
        search_.set_option_value("min_node_level", 0);
        search_.set_option_value("max_solutions", 100000000);
        search_.set_option_value("accept_score", 15);
        search_.set_option_value("max_size", 1000);
        
        for(auto & opt : opts) {
            if(! search_.has_option(opt->name())) { continue; }
            
        }
    }
    
    
private:
    Options options_;
    MotifStateSearch search_;
};

#endif /* defined(__RNAMake__mini_ttr__) */
