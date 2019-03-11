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

base::CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv) {
    
    base::CommandLineOptions cl_opts;
    auto search = MotifStateSearch();
    cl_opts.add_options(search.options());
    cl_opts.add_option("path", String(""), base::OptionType::STRING, false);
    cl_opts.add_option("test_run", false, base::OptionType::BOOL, false);
    cl_opts.add_option("out", String("solutions.top"), base::OptionType::STRING, false);
    cl_opts.add_option("opt_seq", true, base::OptionType::BOOL, false);
    
    cl_opts.parse_command_line(argc, argv);
    return cl_opts;
}


void
MiniTTR::run() {
    //mg_.replace_ideal_helices();
    auto end = mg_.get_available_end(2);
    auto start  = mg_.get_available_end(3);
    //auto end_pos = mg_.get_end(2);
    mg_.increase_level();
    
    auto beads = mg_.beads();
    auto centers = math::Points();
    for(auto const & b : beads) {
        if(b.btype() == BeadType::PHOS) {
            continue;
        }
        centers.push_back(b.center());
    }
    
    search_.setup(start->state(), end->state());
    search_.beads(centers);

    //if(test_run_) {
        auto sol = search_.next();
        if(sol == nullptr) {
            std::cout << "no solutions" << std::endl;
            exit(0);
        }
        auto mt = sol->to_motif_tree();
        mg_.add_motif_tree(mt, 3);
        mg_.add_connection(2, mg_.last_node()->index(), end->name(), "");
        //mg_.write_pdbs();
        if(opt_seq_) {
            optimize_sequence(mg_);
        }
        return;
    //}
    
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
        mg_.add_motif_tree(mt, 3);

        
        mg_.add_connection(2, mg_.last_node()->index(), end->name(), "");
        out << mg_.topology_to_str() << std::endl;
        try{
            mg_.to_pdb("test."+ std::to_string(count) + ".pdb", 1);
        }
        catch(...) { }
        mg_.remove_level(1);

    }
    
    
}

void
MiniTTR::optimize_sequence(MotifGraph & org_mg) {
    auto str = org_mg.topology_to_str();
    auto mg = std::make_shared<MotifGraph>(str, MotifGraphStringType::OLD);
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

    //auto r = optimizer_.optimize(mg, free_end_node, tetraloop_node, 1, 2);
    //std::cout << r->score << std::endl;
    //r->motif_tree->write_pdbs();
    //r->motif_tree->to_pdb("final.pdb");
    
    
}

void
MiniTTRPathFollow::run() {
    
    auto lines =base::get_lines_from_file(options_.get_string("path"));
    auto path_points = math::vectors_from_str(lines[0]);
    auto scorer = std::make_shared<MSS_PathFollow>(path_points);
    
    auto end = mg_.get_available_end(0, "A222-A251");
    auto start  = mg_.get_available_end(mg_.last_node()->index());
    auto end_pos = mg_.last_node()->index();
    
    auto beads = mg_.beads();
    auto centers = math::Points();
    for(auto const & b : beads) {
        if(b.btype() == BeadType::PHOS) {
            continue;
        }
        centers.push_back(b.center());
    }

    search_.setup(start->state(), end->state());
    search_.scorer(scorer);
    search_.beads(centers);
  
    /*auto selector = std::make_shared<MotifStateSelector>(MSS_RoundRobin());
    selector->add("ideal_helices");
    selector->add("unique_twoway");
    search_.selector(selector);
    */

    if(test_run_) {
        search_.set_option_value("verbose", true);
        auto sol = search_.next();
        if(sol == nullptr) {
            std::cout << "no solutions" << std::endl;
            exit(0);
        }
        auto mt = sol->to_motif_tree();
        mg_.add_motif_tree(mt, 1, "A1-A6");
        mg_.write_pdbs();
        mg_.to_pdb("r.pdb", 1);
        exit(0);
        auto start_2 = mg_.get_available_end(mg_.last_node()->index());
        
        auto search_2 = MotifStateSearch();
        search_2.set_option_value("accept_score", 10.0f);
        
        beads = mg_.beads();
        centers = math::Points();
        for(auto const & b : beads) {
            centers.push_back(b.center());
        }
        
        search_2.beads(centers);
        search_2.setup(start_2->state(), end->state());
        auto sol_2 = search_2.next();
        auto mt2 = sol_2->to_motif_tree();
        mg_.add_motif_tree(mt2, mg_.last_node()->index());
        mg_.add_connection(0, mg_.last_node()->index(), end->name(), "");

        mg_.to_pdb("full_path.pdb", 1);
        mg_.write_pdbs();
    }
    
    
}



int main(int argc, const char * argv[]) {
    auto base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/mini_ttr/resources/";
    RM::instance().add_motif(base_path+"GAAA_tetraloop");

    auto options = parse_command_line(argc, argv);
    auto app = std::make_shared<MiniTTR>();
    
    if(options.is_filled("path")) { app.reset(new MiniTTRPathFollow()); }
    
    app->setup(options);
    app->run();
    
    exit(0);
    
    
    return 0;
}























