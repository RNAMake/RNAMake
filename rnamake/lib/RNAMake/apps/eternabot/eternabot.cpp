//
//  eternabot.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "eternabot.h"

#include "base/cl_option.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "eternabot/sequence_designer.h"

CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv) {
    
    CommandLineOptions cl_opts;
    cl_opts.add_option("seq", String("NNNNAAAANNNN"), OptionType::STRING, true);
    cl_opts.add_option("ss",  String("((((....))))"), OptionType::STRING, true);
    cl_opts.add_option("steps", 100, OptionType::INT, false);

    cl_opts.parse_command_line(argc, argv);
    return cl_opts;
}

int main(int argc, const char * argv[]) {
    auto cl_opts = parse_command_line(argc, argv);
    
    auto designer = eternabot::SequenceDesigner();
    
    if(cl_opts.is_filled("steps")) {
        designer.set_option_value("steps", cl_opts.get_int("steps"));
    }
    
    designer.setup();
    auto parser = sstruct::SecondaryStructureParser();
    auto p = parser.parse_to_pose(cl_opts.get_string("seq"),
                                  cl_opts.get_string("ss"));
    auto results = designer.design(p);
    std::cout << results[0]->score << " " << results[0]->sequence << std::endl;
    
}
