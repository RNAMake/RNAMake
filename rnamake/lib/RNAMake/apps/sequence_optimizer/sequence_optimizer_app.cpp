//
//  sequence_optimizer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "util/file_io.h"
#include "sequence_optimizer_app.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "motif_data_structures/motif_graph.h"

void
SequenceOptimizerApp::setup_options() {
    add_option("f", String(""), OptionType::STRING, true);
    add_option("break_point", String(""), OptionType::STRING, true);
    add_option("connections", String(""), OptionType::STRING, false);
    add_option("o", String("sequences.dat"), OptionType::STRING, false);
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
    
    auto file_path = get_string_option("f");
    auto lines = get_lines_from_file(file_path);
    auto connections = std::vector<Ints>();
    auto spl = split_str_by_delimiter(get_string_option("connections"), ",");
    for(auto const & s : spl) {
        auto bounds = split_str_by_delimiter(s, ":");
        connections.push_back(Ints{std::stoi(bounds[0]), std::stoi(bounds[1])});
    }
    
    std::ofstream out;
    out.open(get_string_option("o"));
    
    spl = split_str_by_delimiter(get_string_option("break_point"), ":");
    auto break_point_num = std::stoi(spl[0]);
    auto break_point_name = spl[1];
    int free_ends = 0;
    int free_end_node = 0;
    int i = 0;
    for(auto const & l : lines) {
        auto mg = std::make_shared<MotifGraph>(l, MotifGraphStringType::TOP);
        std::cout << i << std::endl;
        mg->set_option_value("sterics", false);
        for(auto const & c : connections) {
            int start = c[0];
            int end = c[1];
            if(end == -1) { end = mg->last_node()->index(); }
            mg->add_connection(start, end, "", "");
           
        }
        
        mg->replace_ideal_helices();
        
        auto n = mg->get_node(break_point_num);
        auto ei = n->data()->end_index(break_point_name);

        auto c = n->connections()[ei];
        auto partner = c->partner(n->index());
        auto uuid_1 = n->data()->id();
        auto uuid_2 = partner->data()->id();
        auto pei = c->end_index(partner->index());
        
        auto results = optimizer_.get_optimized_sequences_2(mg, uuid_1, uuid_2, ei, pei);
        for(auto const & r : results) {
            out << i << " " << r->close_distance << " " << r->eternabot_score << " " << r->sequence << std::endl;
        }
        i++;
        
    }
    
    out.close();
    
    
}

int main(int argc, const char * argv[]) {
    
    //load TTR
    auto base_dir = String("/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/apps/simulate_tectos");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GAAA_tetraloop");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GGAA_tetraloop");

    
    auto app = SequenceOptimizerApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;

}