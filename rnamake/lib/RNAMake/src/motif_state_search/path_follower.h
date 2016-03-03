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

#include "base/option.h"
#include "base/cl_option.h"

#include "util/settings.h"
#include "motif/motif.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_state_search/motif_state_search.h"

class PathFollower {
public:
    PathFollower() { setup_options(); }
    
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
        
        search_.set_option_value("max_node_level", 8);
        search_.set_option_value("min_node_level", 0);
        search_.set_option_value("max_solutions", 100000000);
        search_.set_option_value("accept_score", -10000.1f);
        search_.set_option_value("max_size", 10000);
        search_.set_option_value("max_steps", 10000.0f);
        search_.set_option_value("verbose", false);
        
        auto beads = mg_->beads();
        auto centers = Points();
        for(auto const & b : beads) {
            if(b.btype() == BeadType::PHOS) {
                continue;
            }
            centers.push_back(b.center());
        }

        auto scorer = std::make_shared<MSS_PathFollow>(path);
        start_state_ = mg_->get_node(ni)->data()->get_basepair(end_name)[0]->state();
        search_.setup(start_state_, start_state_);
        search_.scorer(scorer);
        search_.beads(centers);
    }
    
    void
    setup(
        Points const & path) {
        path_ = path;
        
        //set default search values
        search_.set_option_value("max_node_level", 8);
        search_.set_option_value("min_node_level", 0);
        search_.set_option_value("max_solutions", 100000000);
        search_.set_option_value("accept_score", -10000.1f);
        search_.set_option_value("max_size", 10000);
        search_.set_option_value("max_steps", 10000.0f);
        search_.set_option_value("verbose", false);

        auto scorer = std::make_shared<MSS_PathFollow>(path);
        auto f_path = motif_dirs() + "ref.motif";

        start_ = file_to_motif(f_path)->ends()[0];
        search_.setup(start_->state(), start_->state());
        search_.scorer(scorer);
        
    }
    
    void
    setup(
        Points const & path,
        BasepairStateOP const & start,
        Points const & centers) {
        
        path_ = path;
        start_state_ = start;
        
        //set default search values
        search_.set_option_value("max_node_level", 8);
        search_.set_option_value("min_node_level", 0);
        search_.set_option_value("max_solutions", 100000000);
        search_.set_option_value("accept_score", -10000.1f);
        search_.set_option_value("max_size", 10000);
        search_.set_option_value("max_steps", 10000.0f);
        search_.set_option_value("verbose", false);
        
        auto scorer = std::make_shared<MSS_PathFollow>(path);
        search_.setup(start_state_, start_state_);
        search_.scorer(scorer);
        search_.beads(centers);

    }
    
    void
    set_cmd_options(
        CommandLineOptions const &);
    
    MotifTreeOP 
    next();
    
    MotifStateSearchSolutionOPs
    solutions();
    
    
public: //option wrappers
    
    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }
    
    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }
    
    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }
    
    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }
    
    inline
    bool
    has_option(String const & name) { return options_.has_option(name); }
    
    template<typename T>
    void
    set_option_value(
        String const & name,
        T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }
    
private:
    void
    setup_options();
    
    void
    update_var_options();
    
        
private:
    MotifGraphOP mg_;
    Points path_;
    int ni_;
    BasepairOP start_;
    BasepairStateOP start_state_;
    MotifStateSearch search_;
    Options options_;

    
};

#endif /* defined(__RNAMake__path_follower__) */
