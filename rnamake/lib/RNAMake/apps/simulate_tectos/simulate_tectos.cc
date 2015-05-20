//
//  simulate_tectos.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "simulate_tectos.h"
#include "base/cl_option.h"

Options
parse_command_line(
    int argc,
    const char ** argv) {
    
    CL_Options cl_opts;
    cl_opts.add_option("fseq", "", STRING_TYPE, "CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG", false);
    cl_opts.add_option("fss" , "", STRING_TYPE, "((((((....((((((((((((....))))))))))))....))))))", false);
    cl_opts.add_option("cseq", "", STRING_TYPE, "CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG", false);
    cl_opts.add_option("css" , "", STRING_TYPE, "(((((((..((((((((((((....))))))))))))...)))))))", false);
    cl_opts.add_option("s", "steps", FLOAT_TYPE, "1000000", false);
    
    return cl_opts.parse_command_line(argc, argv);
    
}

int main(int argc, const char * argv[]) {
    try {
        Options opts = parse_command_line(argc, argv);
    } catch(std::runtime_error e) {
        std::cerr << "caught runtime exception: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    return 0;
}