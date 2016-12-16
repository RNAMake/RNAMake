//
//  sequence_optimizer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//


#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "motif_data_structures/motif_graph.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"
#include "sequence_optimizer_app.hpp"


void
SequenceOptimizerApp::setup_options() {
    add_option("mg", String(""), OptionType::STRING);
    
    add_cl_options(optimizer_.options(), "optimizer");
}

void
SequenceOptimizerApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, optimizer_.options(), "optimizer");
}

void
SequenceOptimizerApp::run() {
    
    auto mg = MotifGraphOP();
    
    auto lines = get_lines_from_file(get_string_option("mg"));
    for(auto l : lines) {
        mg = std::make_shared<MotifGraph>(l, MotifGraphStringType::MG);
        /*try {
            mg = std::make_shared<MotifGraph>(l, MotifGraphStringType::MG);
        }
        catch(...) { break; }*/
        
        auto c = GraphtoTree();
        auto d_mt = c.convert(mg, nullptr, -1, nullptr);
        d_mt->set_option_value("sterics", false);
        
        d_mt->write_pdbs();
        
        exit(0);

    }
    
    
    
}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);
    
    //load tectos
    auto tecto_dir = String(base_dir()+"/rnamake/lib/RNAMake/apps/simulate_tectos");
    RM::instance().add_motif(tecto_dir+"/resources/GAAA_tetraloop");
    RM::instance().add_motif(tecto_dir+"/resources/GGAA_tetraloop");

    
    auto app = SequenceOptimizerApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;

}
