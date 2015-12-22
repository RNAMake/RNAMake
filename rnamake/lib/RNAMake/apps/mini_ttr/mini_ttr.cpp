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
    
    auto lines = get_lines_from_file("all_points.str");
    auto path_points = vectors_from_str(lines[0]);
 
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
    search.option("max_node_level", 20);
    search.option("min_node_level", 6);
    search.option("max_solutions", 100000000);
    search.option("accept_score", 20.0f);
    search.option("max_size", 500);
    //search.option("min_ss_score", 0.0f);
    search.setup(start->state(), end->state());
    search.beads(centers);
    search.path(path_points);
    
    std::ofstream out;
    out.open("solutions.top", std::ofstream::out);
    
    mg.option("sterics", 0);
    mg.increase_level();
    int count = 0;
    while(! search.finished()) {
        auto sol = search.next();
        if(sol == nullptr) {
            std::cout << "no solutions" << std::endl;
            exit(0);
        }
        count++;
        if(count % 100 == 0) {
            std::cout << count << std::endl;
        }
        
        auto mt = sol->to_motif_tree();
        mg.add_motif_tree(mt, 0, "A222-A251");
        mg.add_connection(end_pos, mg.last_node()->index(), end->name());
        mg.replace_ideal_helices();
        mg.write_pdbs();
        exit(0);
        out << mg.topology_to_str() << std::endl;
        
        mg.remove_level(1);

        if(count % 1000 == 0) {
            mg = MotifGraph();
            mg.add_motif("GAAA_tetraloop", "A229-A245");
            mg.add_motif("HELIX.IDEAL.6", -1, "A149-A154");
            mg.replace_ideal_helices();
            mg.increase_level();

        }
        
    }
    
    out.close();
    
    
    //sol->to_pdb("test.pdb", 1);
    
    
    return 0;
}

