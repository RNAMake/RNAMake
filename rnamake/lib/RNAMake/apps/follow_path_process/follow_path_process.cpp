//
//  follow_path_process.c
//  RNAMake
//
//  Created by Joseph Yesselman on 3/10/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "follow_path_process.h"
#include "base/cl_option.h"
#include "util/file_io.h"
#include "util/basic_io.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "motif_data_structures/motif_graph.h"
#include "sequence_optimizer/sequence_optimizer.h"


CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv) {
    
    CommandLineOptions cl_opts;
    cl_opts.add_option("f", String(""), OptionType::STRING, true);
    cl_opts.parse_command_line(argc, argv);
    
    return cl_opts;
    
}


int main(int argc, const char * argv[]) {
    //auto cmd_opts = parse_command_line(argc, argv);
    
    //load TTR
    auto base_dir = String("/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/apps/mini_ttr");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GAAA_tetraloop");
    
    auto lines = get_lines_from_file("solutions.top");
    auto mt = std::make_shared<MotifTree>(lines[0]);
    auto mg = std::make_shared<MotifGraph>();
    mg->set_option_value("sterics", false);
    mg->add_motif_tree(mt);
    mg->add_connection(3, mg->last_node()->index(), "", "");
    
    std::cout << "START!" << std::endl;
    //std::cout << mg->secondary_structure()->sequence() << std::endl;
    std::cout << mg->secondary_structure()->dot_bracket() << std::endl;
    exit(0);
    
    mg->replace_ideal_helices();
    
    int tetraloop_node = -1, free_end_node = -1;
    int free_ends = 0;
    
    for(auto const & n : *mg) {
        if(n->data()->name() == "GAAA_tetraloop") {
            tetraloop_node = n->index();
            continue;
        }
        
        free_ends = 0;
        for(auto const & c : n->connections()) {
            if(c == nullptr) { free_ends++; }
        }
        if(free_ends > 0) {
            free_end_node = n->index();
        }
        
    }
    
    
    auto seq_opt = SequenceOptimizer();
    seq_opt.optimize(mg, free_end_node, tetraloop_node, 0, 2);
}



















