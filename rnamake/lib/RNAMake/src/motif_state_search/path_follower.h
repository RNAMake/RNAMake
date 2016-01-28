//
//  path_follower.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/18/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__path_follower__
#define __RNAMake__path_follower__

#include <stdio.h>


#include "motif_data_structures/motif_graph.h"
#include "motif_state_search/motif_state_search.h"

class PathFollower {
public:
    PathFollower() {}
    
    ~PathFollower() {}
        
public:
    
    void
    setup(
        Points const & path,
        MotifGraphOP const & mg,
        int ni,
        String const & end_name) {
        path_ = path;
        mg_ = mg;
        ni_ = ni;
        
        search_.set_option_value("max_node_level", 400);
        search_.set_option_value("min_node_level", 0);
        search_.set_option_value("max_solutions", 100000000);
        search_.set_option_value("accept_score", 0.1f);
        search_.set_option_value("max_size", 10000);
        search_.set_option_value("max_steps", 100000.0f);
        search_.set_option_value("verbose", true);
        
        auto beads = mg_->beads();
        auto centers = Points();
        for(auto const & b : beads) {
            if(b.btype() == BeadType::PHOS) {
                continue;
            }
            centers.push_back(b.center());
        }

        auto scorer = std::make_shared<MSS_PathFollow>(path);
        start_ = mg_->get_node(ni)->data()->get_basepair(end_name)[0];
        search_.setup(start_->state(), start_->state());
        search_.scorer(scorer);
        search_.beads(centers);
    }
    
    MotifTreeOP 
    next();
    
        
private:
    MotifGraphOP mg_;
    Points path_;
    int ni_;
    BasepairOP start_;
    MotifStateSearch search_;
        
};

#endif /* defined(__RNAMake__path_follower__) */
