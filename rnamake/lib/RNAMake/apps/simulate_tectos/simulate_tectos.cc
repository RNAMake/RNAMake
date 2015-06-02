//
//  simulate_tectos.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/cl_option.h"
#include "util/settings.h"
#include "secondary_structure/ss_tree.h"
#include "motif/motif_tree_topology.h"
#include "resources/library_manager.h"
#include "simulate_tectos.h"

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

SimulateTectos::SimulateTectos(
    String const & fseq,
    String const & fss,
    String const & cseq,
    String const & css) {
    
    int f_diff = ((int)fseq.length()-10)-10;
    String sub_fseq = fseq.substr(10, f_diff);
    String sub_fss  = fss.substr(10, f_diff);
    
    //std::cout << sub_fseq << std::endl;
    //std::cout << sub_fss << std::endl;
    
    SS_Tree flow_tree(sub_fss, sub_fseq);
    MotifTreeTopology mtt(flow_tree);
    //SS_Tree chip_tree(css, cseq);
    
    

}


int main(int argc, const char * argv[]) {
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    LibraryManager::getInstance().add_motif(base_path+"GAAA_tetraloop");
    LibraryManager::getInstance().add_motif(base_path+"GGAA_tetraloop");

    try {
        
        Options opts = parse_command_line(argc, argv);
        SimulateTectos st(opts.option<String>("fseq"),
                          opts.option<String>("fss"),
                          opts.option<String>("cseq"),
                          opts.option<String>("css"));
                                              
    } catch(std::runtime_error e) {
        std::cerr << "caught runtime exception: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    return 0;
}