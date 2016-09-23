//
//  sequence_optimizer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//


#include "sequence_optimizer_app.hpp"
#include "util/steric_lookup.hpp"
#include "base/file_io.h"
#include "base/settings.h"
#include "motif/motif_factory.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "motif_data_structures/motif_graph.h"
#include "eternabot/sequence_designer.h"

void
SequenceOptimizerApp::setup_options() {
    add_option("f", String(""), OptionType::STRING, true);
    add_option("break_point", String(""), OptionType::STRING, false);
    add_option("connections", String(""), OptionType::STRING, false);
    add_option("o", String("sequences.dat"), OptionType::STRING, false);
    
    
    add_option("mg", String(""), OptionType::STRING);
    
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
    
    if(get_string_option("mg") != "") {
        run_from_mg();
    }
    
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
        
        auto results = optimizer_.get_optimized_sequences(mg, uuid_1, uuid_2, ei, pei);
        for(auto const & r : results) {
            out << i << " " << r->close_distance << " " << r->eternabot_score << " " << r->sequence << std::endl;
        }
        i++;
        
    }
    
    out.close();
    
    
}

void
SequenceOptimizerApp::run_from_mg() {
    
    auto file_path = get_string_option("mg");
    auto lines = get_lines_from_file(file_path);
    
    auto mg = std::make_shared<MotifGraph>(lines[0], MotifGraphStringType::MG);
    mg->write_pdbs("ribosome");
    auto spl = split_str_by_delimiter(lines[1], " ");
    auto start_name = spl[0];
    auto start_pos = std::stoi(spl[1]);
    spl = split_str_by_delimiter(lines[2], " ");
    auto end_name = spl[0];
    auto end_pos = std::stoi(spl[1]);
    
    auto start = mg->get_node(start_pos)->data()->get_basepair(start_name)[0];
    auto end = mg->get_node(end_pos)->data()->get_basepair(end_name)[0];

    auto mf = MotifFactory();
    start->bp_type("cW-W");
    
    auto res = mg->get_node(start_pos)->data()->get_residue(start->res1()->num()+1, "A", "");
    
    auto bp_2 = mg->get_node(start_pos)->data()->get_basepair(res->uuid())[0];
    bp_2->bp_type("cW-W");
    
    auto m = mf.motif_from_bps({bp_2, start});
    
    //m->block_end_add(0);
    
    auto mg_build = std::make_shared<MotifGraph>();
    mg_build->add_motif(m);
    mg_build->increase_level();
    //mg_build->add_motif(get_motif_from_resource_manager("HELIX.IDEAL.2"));
    //mg_build->write_pdbs();
    
    file_path = get_string_option("f");
    auto sol_lines = get_lines_from_file(file_path);
    
    auto cols = split_str_by_delimiter(sol_lines[0], "\t");
    std::cout << cols[1] << std::endl;
    
    auto beads = Points();
    for(auto & n : *mg) {
        //n->data()->get_beads(n->data()->ends());
        for(auto const & b : n->data()->beads()) {
            if(b.btype() == BeadType::PHOS) { continue; }
            beads.push_back(b.center());
            
        }
        
    }
    
    std::ofstream out;
    out.open(get_string_option("o"));
    
    auto sl = StericLookup();
    sl.add_points(beads);
    
    auto end_state_1 = end->state();
    auto end_state_2 = BasepairStateOP();
    auto end_state_2_flip = BasepairStateOP();
    int clash = 0;
    float dist = 0;
    
    int i = -1;
    for(auto const & l : sol_lines) {
        i++;
        if(i == 0) { continue; }
        auto sol_spl = split_str_by_delimiter(l, "\t");
        auto mt = std::make_shared<MotifTree>(sol_spl[4]);
        mg_build->add_motif_tree(mt, 0, start_name);
        mg_build->replace_ideal_helices();

        auto designer = eternabot::SequenceDesigner();
        designer.set_option_value("steps", 1000);
        designer.setup();
        auto dss = mg_build->designable_secondary_structure();
        
        auto results = designer.design(dss);
        int j = 0;
        for(auto const & r : results) {
            j++;
            //std::cout << r->sequence << " " << r->score << std::endl;
            dss->replace_sequence(r->sequence);
            mg_build->replace_helical_sequence(dss);
            auto new_points = Points();
            for(auto const & n : *mg_build) {
                for(auto const & b : n->data()->beads()) {
                    if(b.btype() == BeadType::PHOS) { continue; }
                    new_points.push_back(b.center());
                }
            }
            
            clash = sl.clash(new_points);
            end_state_2 = mg_build->last_node()->data()->ends()[1]->state();
            end_state_2->flip();
            end_state_2_flip = std::make_shared<BasepairState>(end_state_2->copy());
            end_state_2->flip();
            
            dist = new_score_function_new(end_state_1, end_state_2, end_state_2_flip);

            out << i << " " << j << " " << dist << " " << clash << " " << r->sequence << " " << r->score << std::endl;
            std::cout << i << " " << j << " " << dist << " " << clash << " " << r->sequence << " " << r->score << std::endl;
            
            
            //rebuild to avoid weird behavior
            mg_build = std::make_shared<MotifGraph>();
            mg_build->add_motif(m);
            mg_build->add_motif_tree(mt, 0, start_name);
            mg_build->replace_ideal_helices();
            dss = mg_build->designable_secondary_structure();
            
            //exit(0);
            
            if(j > 100) {
                break;
            }
            
        }
        
        mg_build = std::make_shared<MotifGraph>();
        mg_build->add_motif(m);
        
    }
    
    out.close();
  
    
    
    //design_mg
    
    
    exit(0);
}

int main(int argc, const char * argv[]) {
    
    //load tectos
    auto tecto_dir = String(base_dir()+"/rnamake/lib/RNAMake/apps/simulate_tectos");
    RM::instance().add_motif(tecto_dir+"/resources/GAAA_tetraloop");
    RM::instance().add_motif(tecto_dir+"/resources/GGAA_tetraloop");

    
    auto app = SequenceOptimizerApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;

}
