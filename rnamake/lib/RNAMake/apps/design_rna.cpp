//
//  design_rna.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/26/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "design_rna.hpp"

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"



void
DesignRNAApp::setup_options() {
    // start from pdb
    add_option("pdb", String(""), OptionType::STRING, false);
    add_option("start_bp", String(""), OptionType::STRING, false);
    add_option("end_bp", String(""), OptionType::STRING, false);

    add_cl_options(search_.options(), "search");
}

void
DesignRNAApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, search_.options(), "search");
    search_.update_var_options();

}


void
DesignRNAApp::_setup_from_pdb() {
    auto struc = RM::instance().get_structure(get_string_option("pdb"), "scaffold");
    
    std::cout << struc->name() << std::endl;
}


void
DesignRNAApp::run() {
    
    if(get_string_option("pdb") != "") { _setup_from_pdb(); }
    
}


int main(int argc, const char * argv[]) {
    std::set_terminate(print_backtrace);
    
    auto app = DesignRNAApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    return 0;
    
}