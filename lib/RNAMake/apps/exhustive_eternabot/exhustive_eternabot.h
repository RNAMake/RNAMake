//
//  exhustive_eternabot.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/5/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__exhustive_eternabot__
#define __RNAMake__exhustive_eternabot__

#include <stdio.h>

#include "base/option.h"
#include "base/cl_option.h"
#include "secondary_structure/pose.h"

base::CommandLineOptions
parse_command_line(int, const char **);

class ExhustiveEternabot {
public:
    ExhustiveEternabot():
    enumerated_bps_(secondary_structure::BasepairOPs()),
    out_name_("4bps.out")
    {}
    
    ~ExhustiveEternabot() {}
    
public:
    
    void
    setup(base::CommandLineOptions const &);
    
    void
    run();

    
    
private:
    std::vector<Strings> pairs_;
    secondary_structure::PoseOP p_;
    secondary_structure::BasepairOPs enumerated_bps_;
    String out_name_;
    
};



#endif /* defined(__RNAMake__exhustive_eternabot__) */
