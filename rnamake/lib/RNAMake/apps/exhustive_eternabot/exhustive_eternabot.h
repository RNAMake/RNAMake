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

CommandLineOptions
parse_command_line(int, const char **);

class ExhustiveEternabot {
public:
    ExhustiveEternabot():
    enumerated_bps_(sstruct::BasepairOPs())
    {}
    
    ~ExhustiveEternabot() {}
    
public:
    
    void
    setup(CommandLineOptions const &);
    
    void
    run();
    
    
private:
    std::vector<Strings> pairs_;
    sstruct::PoseOP p_;
    sstruct::BasepairOPs enumerated_bps_;
    
};



#endif /* defined(__RNAMake__exhustive_eternabot__) */
