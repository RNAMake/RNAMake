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
#include "util/settings.h"
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
    
    //load TTR
    auto ttr_dir = String(base_dir()+"/rnamake/lib/RNAMake/apps/mini_ttr");
    ResourceManager::getInstance().add_motif(ttr_dir+"/resources/GAAA_tetraloop");
    
    std::ofstream out;
    out.open("finished.out");
    
    int i = 0;
    auto lines = get_lines_from_file("solutions.top");
    for(auto const & l : lines) {
        
        std::cout << i << std::endl;
        auto mt = std::make_shared<MotifTree>(l);
        auto mg = std::make_shared<MotifGraph>();
        mg->set_option_value("sterics", false);
        mg->add_motif_tree(mt);
        mg->add_connection(3, mg->last_node()->index(), "", "");
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
        auto sols = seq_opt.get_optimized_sequences(mg, free_end_node, tetraloop_node, 0, 2);
        for(auto const & sol : sols) {
            out << i << " " << sol->close_distance << " " << sol->eternabot_score << " " << sol->sequence << std::endl;
        }
        i++;

    }
    
    out.close();
    
}



















