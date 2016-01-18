//
//  motif_topology_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_topology.h"
#include "util/settings.h"
#include "util/file_io.h"
#include "resources/resource_manager.h"

#include "motif_topology_unittests.h"
#include "build/build_motif_graph.h"


namespace unittests {
namespace motif_structures {

void
MotifTopologyUnittest::test_to_tree() {
    auto builder = BuildMotifGraph();
    auto mg = builder.build(3);
    auto mt = graph_to_tree(mg);
    if(mg->size() != mt->size()) {
        throw UnittestException("did not produce a motif tree of the right size");
    }
    auto mg_end_d = mg->get_node(2)->data()->ends()[1]->d();
    auto mt_end_d = mt->get_node(2)->data()->ends()[1]->d();
    
    auto dist = mg_end_d.distance(mt_end_d);
    if(dist > 0.1) {
        throw UnittestException("graph and tree did not end in the same place");
    }
}
    
void
MotifTopologyUnittest::test_to_tree_complex() {
    auto path = unittest_resource_dir() + "/motif_data_structures/graph_to_tree.top";
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
        std::cout << n->data()->name() << std::endl;
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
    
    auto mt = graph_to_tree(mg, mg->get_node(free_end_node),
                            mg->get_node(tetraloop_node)->data()->ends()[2]);
    
}

int
MotifTopologyUnittest::run() {
    test_to_tree();
    test_to_tree_complex();
    return 1;
}

int
MotifTopologyUnittest::run_all() {
    return 1;
}

}
}














