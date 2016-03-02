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


CommandLineOptions
parse_command_line(int, const char **);

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
    
};


#endif /* defined(__RNAMake__follow_path__) */
