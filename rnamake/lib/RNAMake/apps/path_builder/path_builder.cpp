//
//  path_builder.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 3/23/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "path_builder.hpp"

#include "base/cl_option.h"
#include "util/file_io.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_state_search/motif_state_search.h"
#include "eternabot/sequence_designer.h"

CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv) {
    
    CommandLineOptions cl_opts;
    cl_opts.add_option("mg", String(""), OptionType::STRING, true);
    cl_opts.parse_command_line(argc, argv);
    
    return cl_opts;
    
}


int main(int argc, const char * argv[]) {
    auto cmd_opts = parse_command_line(argc, argv);
    auto base_dir = String("/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/apps/mini_ttr");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GAAA_tetraloop");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GGAA_tetraloop");

    
    auto lines = get_lines_from_file(cmd_opts.get_string("mg"));
    if(lines.size() < 3) {
        throw std::runtime_error("motif graph file did not contain three lines, first line is" \
                                 "topology followed by starting basepair and end basepair");
    }
    
    auto mg = std::make_shared<MotifGraph>(lines[0]);
    auto spl1 = split_str_by_delimiter(lines[1], " ");
    auto spl2 = split_str_by_delimiter(lines[2], " ");

    auto start = mg->get_node(std::stoi(spl1[0]))->data()->get_basepair(spl1[1])[0];
    auto end   = mg->get_node(std::stoi(spl2[0]))->data()->get_basepair(spl2[1])[0];

    auto beads = mg->beads();
    auto centers = Points();
    for(auto const & b : beads) {
        if(b.btype() != BeadType::PHOS) {
            centers.push_back(b.center());
        }
    }
    
    auto search = MotifStateSearch();
    search.setup(start->state(), end->state());
    search.beads(centers);
    mg->increase_level();
    mg->set_option_value("sterics", false);
    mg->write_pdbs();
    auto designer = eternabot::SequenceDesigner();
    
    std::ofstream out;
    out.open("solutions.top");
    
    int i = 1;
    while(!search.finished() ) {
        auto sol = search.next();
        auto mt_sol = sol->to_motif_tree();
        mg->add_motif_tree(mt_sol, std::stoi(spl1[0]), spl1[1]);
        //mg->add_connection(std::stoi(spl1[0]), mg->last_node()->index(), spl1[1]);
        //mg->replace_ideal_helices();
        out << mg->topology_to_str() << std::endl;
        
        mg->remove_level(1);
        
        i += 1;
        
    }
    
    out.close();
    
    
    return 0;

}