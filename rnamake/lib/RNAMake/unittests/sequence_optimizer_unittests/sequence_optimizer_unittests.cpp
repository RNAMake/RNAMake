//
//  sequence_optimizer_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/11/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_topology.h"
#include "util/settings.h"
#include "util/file_io.h"
#include "resources/resource_manager.h"
#include "sequence_optimizer/sequence_optimizer.h"

#include "sequence_optimizer_unittests.h"

namespace unittests {
namespace sequence_optimizer {
    
int
SequenceOptimizerUnittest::test_creation() {
    auto optimizer = SequenceOptimizer();
    return 0;
}

int
SequenceOptimizerUnittest::test_optimize() {
    /*auto path = unittest_resource_dir() + "/motif_data_structures/graph_to_tree.top";
    auto lines = get_lines_from_file(path);
    
    auto base_dir = String("/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/apps/mini_ttr");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GAAA_tetraloop");
    
    auto mg = std::make_shared<MotifGraph>(lines[0]);
    
    mg->replace_ideal_helices();
    
    
    int free_end_node = 0;
    int free_ends = 0;
    int tetraloop_node = 0;
    for(auto const & n : *mg) {
        free_ends = 0;
        for(auto const & c : n->connections()) {
            if(c == nullptr) { free_ends += 1; }
        }
        if(free_ends) {
            free_end_node = n->index();
        }
        
        if(n->data()->name() == "GAAA_tetraloop") {
            tetraloop_node = n->index();
        }
    }
    
    auto optimizer = SequenceOptimizer();
    optimizer.optimize(mg, free_end_node, tetraloop_node, 1, 2);*/

    
    return 0;
    
}
    
int
SequenceOptimizerUnittest::run() {
    test_creation();
    test_optimize();
    return 0;
}

int
SequenceOptimizerUnittest::run_all() {
    return 0;
}
    
}
}