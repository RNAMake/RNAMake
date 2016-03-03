//
//  follow_path.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/21/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__follow_path__
#define __RNAMake__follow_path__

#include <stdio.h>

#include "base/option.h"
#include "base/cl_option.h"
#include "motif_state_search/path_follower.h"
#include "motif_data_structures/motif_state_tree.h"
#include "motif_state_search/motif_state_search.h"


CommandLineOptions
parse_command_line(int, const char **);

struct PathBuilderNode {
    MotifStateTreeOP mst;
    float path_score, ss_score, diversity_score;
    
    inline
    PathBuilderNode(
        MotifStateTreeOP const & nmst,
        float npath_score,
        float nss_score,
        float ndiversity_score):
    mst(nmst),
    path_score(npath_score),
    ss_score(nss_score),
    diversity_score(ndiversity_score) {
        mst->set_option_value("sterics", false);
    }
    
    inline
    PathBuilderNode(
        PathBuilderNode const & n):
    mst(std::make_shared<MotifStateTree>(*n.mst)),
    path_score(n.path_score),
    ss_score(n.ss_score),
    diversity_score(n.diversity_score)
    {}
    
    void
    add_solution(
        MotifStateSearchSolutionOP const & sol) {
        
        auto sol_mst = sol->to_mst();
        mst->add_mst(sol_mst, -1, 1, "", true);
        add_score(sol_mst, sol->score());
    }
    
    void
    add_score(
        MotifStateTreeOP const & mst,
        int const n_path_score) {
        
        path_score += n_path_score;
        for(auto const & n : *mst) {
            ss_score += n->data()->ref_state->score();
        }
        
    }
};

//to allow sorting of nodes
struct PathBuilderNode_LessThanKey {
    inline
    bool
    operator() (
                PathBuilderNode const & n1,
                PathBuilderNode const & n2) {
        return n1.path_score < n2.path_score;
    }
};

typedef std::vector<PathBuilderNode> PathBuilderNodes;

class PathBuilder {
public:
    PathBuilder() { setup_options(); }
    
    ~PathBuilder() {}
    
public:
    
    void
    setup(
        CommandLineOptions const &);
    
    void
    build();
    
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
    
    std::vector<Points>
    _get_segments(
        Points const &);
    
    std::vector<Points>
    _get_sub_pathes(
        std::vector<Points> const &);

private:
    void
    setup_options();
    
    void
    update_var_options() {}
    
    
private:
    Options options_;
    PathBuilderNodes nodes_;
    
};


#endif /* defined(__RNAMake__follow_path__) */
