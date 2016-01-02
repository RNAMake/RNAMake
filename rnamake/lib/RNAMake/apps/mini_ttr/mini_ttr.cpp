//
//  mini_ttr.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/16/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "mini_ttr.h"
#include "base/cl_option.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_state_search/motif_state_search.h"
#include "motif_state_search/motif_state_search_scorer.h"

CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv) {
    
    CommandLineOptions cl_opts;
    auto search = MotifStateSearch();
    cl_opts.add_options(search.options());
    cl_opts.add_option("path", String(""), OptionType::STRING, false);
    cl_opts.add_option("test_run", false, OptionType::BOOL, false);
    cl_opts.add_option("out", String("solutions.top"), OptionType::STRING, false);
    
    cl_opts.parse_command_line(argc, argv);
    return cl_opts;
}


void
MiniTTR::run() {
    //mg_.replace_ideal_helices();
    auto end = mg_.get_end(2);
    auto start  = mg_.get_end(3);
    //auto end_pos = mg_.get_end(2);
    mg_.increase_level();
    
    auto beads = mg_.beads();
    auto centers = Points();
    for(auto const & b : beads) {
        if(b.btype() == BeadType::PHOS) {
            continue;
        }
        centers.push_back(b.center());
    }
    
    search_.setup(start->state(), end->state());
    search_.beads(centers);

    if(test_run_) {
        auto sol = search_.next();
        if(sol == nullptr) {
            std::cout << "no solutions" << std::endl;
            exit(0);
        }
        auto mt = sol->to_motif_tree();
        mg_.add_motif_tree(mt, 3);
        mg_.write_pdbs();
        return;
    }
    
    std::ofstream out;
    out.open(options_.get_string("out"), std::ofstream::out);
    
    int count = 0;
    while(! search_.finished()) {
        auto sol = search_.next();
        if(sol == nullptr) {
            std::cout << "no more solutions" << std::endl;
            exit(0);
        }
        count++;
        if(count % 100 == 0) {
            std::cout << count << std::endl;
        }
        
        auto mt = sol->to_motif_tree();
        //mg_.add_motif_tree(mt, mg_.last_node()->index());
        mg_.add_motif_tree(mt, 3);

        
        mg_.add_connection(2, mg_.last_node()->index(), end->name());
        //mg_.replace_ideal_helices();
        //mg_.replace_ideal_helices();
        out << mg_.topology_to_str() << std::endl;
        try{
            mg_.to_pdb("test."+ std::to_string(count) + ".pdb", 1);
        }
        catch(...) { }
        mg_.remove_level(1);

    }
    
    
}

void
MiniTTRPathFollow::run() {
    
    auto lines = get_lines_from_file(options_.get_string("path"));
    auto path_points = vectors_from_str(lines[0]);
    auto scorer = std::make_shared<MSS_PathFollow>(path_points);
    
    auto end = mg_.get_end(0, "A222-A251");
    auto start  = mg_.get_end(mg_.last_node()->index());
    auto end_pos = mg_.last_node()->index();
    
    auto beads = mg_.beads();
    auto centers = Points();
    for(auto const & b : beads) {
        if(b.btype() == BeadType::PHOS) {
            continue;
        }
        centers.push_back(b.center());
    }

    search_.setup(start->state(), end->state());
    search_.scorer(scorer);
    search_.beads(centers);
    
    if(test_run_) {
        auto sol = search_.next();
        if(sol == nullptr) {
            std::cout << "no solutions" << std::endl;
            exit(0);
        }
        auto mt = sol->to_motif_tree();
        mg_.add_motif_tree(mt, 1, "A1-A6");
        auto start_2 = mg_.get_end(mg_.last_node()->index());
        
        auto search_2 = MotifStateSearch();
        search_2.setup(start_2->state(), end->state());
        auto sol_2 = search_2.next();
        auto mt2 = sol_2->to_motif_tree();
        mg_.add_motif_tree(mt2, mg_.last_node()->index());
        mg_.write_pdbs();
    }
    
    
}



int main(int argc, const char * argv[]) {
    auto base_path = base_dir() + "/rnamake/lib/RNAMake/apps/mini_ttr/resources/";
    ResourceManager::getInstance().add_motif(base_path+"GAAA_tetraloop");

    auto options = parse_command_line(argc, argv);
    auto app = std::make_shared<MiniTTR>();
    
    if(options.is_filled("path")) { app.reset(new MiniTTRPathFollow()); }
    
    app->setup(options);
    app->run();
    
    exit(0);
    
    
    auto lines = get_lines_from_file("all_points.str");
    auto path_points = vectors_from_str(lines[0]);
 
    auto mg = MotifGraph();
    mg.add_motif("GAAA_tetraloop", "A229-A245");
    mg.add_motif("HELIX.IDEAL.6", -1, "A149-A154");
    //mg.replace_ideal_helices();
    auto end = mg.get_end(0, "A222-A251");
    auto start   = mg.get_end(mg.last_node()->index());
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
    search.set_option_value("max_node_level", 40);
    search.set_option_value("min_node_level", 0);
    search.set_option_value("max_solutions", 100000000);
    search.set_option_value("accept_score", 15);
    search.set_option_value("max_size", 1000);
    //search.set_option_value("sterics", false);
    //search.option("min_ss_score", 0.0f);
    search.setup(start->state(), end->state());
    //search.beads(centers);
    
    auto scorer = std::make_shared<MSS_PathFollow>(path_points);
    search.scorer(scorer);
    //search.path(path_points);
    
    std::ofstream out;
    out.open("solutions.top", std::ofstream::out);
    
    mg.set_option_value("sterics", false);
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
        //mg.add_motif_tree(mt, 1, "A222-A251");
        mg.add_motif_tree(mt, 1, "A1-A6");

        //mg.add_connection(end_pos, mg.last_node()->index(), end->name());
        //mg.replace_ideal_helices();
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
            mg.set_option_value("sterics", false);

        }
        
    }
    
    out.close();
    
    
    //sol->to_pdb("test.pdb", 1);
    
    
    return 0;
}

