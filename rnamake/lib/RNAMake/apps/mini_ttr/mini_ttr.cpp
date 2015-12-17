//
//  mini_ttr.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/16/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "mini_ttr.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_state_search/motif_state_search.h"



int main(int argc, const char * argv[]) {
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/mini_ttr/resources/";
    ResourceManager::getInstance().add_motif(base_path+"GAAA_tetraloop");
    
    auto mg = MotifGraph();
    mg.add_motif("GAAA_tetraloop", "A229-A245");
    mg.add_motif("HELIX.IDEAL.6", -1, "A149-A154");
    mg.replace_ideal_helices();
    auto start = mg.get_end(0, "A222-A251");
    auto end   = mg.get_end(mg.last_node()->index());
    int end_pos = mg.last_node()->index();
    
    auto beads = mg.beads();
    auto centers = Points();
    for(auto const & b : beads) {
        if(b.btype() == BeadType::PHOS) {
            continue;
        }
        centers.push_back(b.center());
    }
    
    auto search = MotifStateSearch();
    search.option("max_node_level", 8);
    search.option("max_solutions", 100000000);
    search.option("accept_score", 10.0f);
    search.setup(start->state(), end->state());
    search.beads(centers);
    
    auto sol = search.next();
    if(sol == nullptr) {
        std::cout << "no solution found" << std::endl;
        exit(0);
    }
    auto mt = sol->to_motif_tree();
    mg.option("sterics", 0);
    mg.add_motif_tree(mt, 0, "A222-A251");
    
    
    exit(0);
    
    int count = 0;
    while(! search.finished()) {
        auto sol = search.next();
        count++;
        if(count % 100 == 0) {
            std::cout << count << std::endl;
        }
        
        //break;
    }
    
    
    
    //sol->to_pdb("test.pdb", 1);
    
    
    return 0;
}

