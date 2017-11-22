//
//  sequence_optimizer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <stdexcept>

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "motif_data_structures/motif_graph.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"
#include "sequence_optimizer_app.hpp"


void
SequenceOptimizerApp::setup_options() {
    add_option("mg", String(""), OptionType::STRING, true);
    add_option("end_1", String(""), OptionType::STRING, true);
    add_option("end_2", String(""), OptionType::STRING, true);
    
    
    add_option("v", false, OptionType::BOOL);
    add_option("out_file", "default.out", OptionType::STRING);
    add_option("score_file", "default.scores", OptionType::STRING);
    add_option("n", 1, OptionType::INT);
    add_option("opt", String("Internal"), OptionType::STRING);
    add_option("pdbs", false, OptionType::BOOL);
    
    add_cl_options(optimizer_.options(), "optimizer");
}

void
SequenceOptimizerApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, optimizer_.options(), "optimizer");
    optimizer_.set_option_value("verbose", get_bool_option("v"));
}

void
SequenceOptimizerApp::run() {
    
    std::ofstream out, sf_out;
    sf_out.open(get_string_option("score_file"));
    sf_out << "opt_num,opt_score,eterna_score,opt_sequence,opt_structure" << std::endl;
    
    out.open(get_string_option("out_file"));

    auto opt_num = 1;
    auto lines = get_lines_from_file(get_string_option("mg"));
    auto mg = std::make_shared<MotifGraph>(lines[0], MotifGraphStringType::MG);
    auto spl = split_str_by_delimiter(get_string_option("end_1"), ",");
    
    if(spl.size() != 2) {
        throw std::runtime_error(
            "incorrect format for end_1 must be graph_pos,end_pos, example 20,1 for the 20th node "
            "and end 1");
    }
    auto ni1 = std::stoi(spl[0]);
    auto ei1 = std::stoi(spl[1]);
    
    spl = split_str_by_delimiter(get_string_option("end_2"), ",");
    if(spl.size() != 2) {
        throw std::runtime_error(
            "incorrect format for end_2 must be graph_pos,end_pos, example 20,1 for the 20th node "
            "and end 1");
    }
    
    auto ni2 = std::stoi(spl[0]);
    auto ei2 = std::stoi(spl[1]);
    auto scorer = SequenceOptimizerScorerOP(nullptr);
    
    if(get_string_option("opt") == "Internal") {
        scorer = std::make_shared<InternalTargetScorer>(ni1, ei1, ni2, ei2);
    }
    
    auto mg_copy = std::make_shared<MotifGraph>(*mg);
    auto write_pdbs = get_bool_option("pdbs");
    
    for(int i = 0; i < get_int_option("n"); i++) {
        auto sols = optimizer_.get_optimized_sequences(mg, scorer);
        for(auto const & s : sols ) {
            mg_copy->replace_helical_sequence(s->sequence);
            sf_out << opt_num << "," << s->dist_score << "," << s->eterna_score << "," << s->sequence;
            sf_out << "," << mg_copy->dot_bracket() << std::endl;
            out << mg_copy->to_str() << std::endl;
            if(write_pdbs) {
                try {
                    mg_copy->to_pdb("design."+std::to_string(i)+".pdb", 1);
                }
                catch(...) { continue; }
            }
            //mg_copy->write_pdbs();
            opt_num += 1;
        
        }
    }
    
    sf_out.close();
    out.close();
    
    
}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);
    
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
