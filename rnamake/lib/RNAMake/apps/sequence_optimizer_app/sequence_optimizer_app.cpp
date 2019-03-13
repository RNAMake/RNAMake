//
//  sequence_optimization.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/2/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <stdexcept>

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_topology.h"
#include "motif_data_structure/motif_graph.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"
#include "sequence_optimizer_app.hpp"


void
SequenceOptimizerApp::setup_options() {
    add_option("mg", String(""), base::OptionType::STRING, true);
    add_option("end_1", String(""), base::OptionType::STRING, false);
    add_option("end_2", String(""), base::OptionType::STRING, false);

    add_option("connections", String(""), base::OptionType::STRING, false);
    
    
    add_option("v", false, base::OptionType::BOOL);
    add_option("out_file", "default.out", base::OptionType::STRING);
    add_option("score_file", "default.scores", base::OptionType::STRING);
    add_option("n", 1, base::OptionType::INT);
    add_option("opt", String("Internal"), base::OptionType::STRING);
    add_option("pdbs", false, base::OptionType::BOOL);
    
    add_cl_options(optimizer_.options(), "optimizer");
}

void
SequenceOptimizerApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    base::Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, optimizer_.options(), "optimizer");
    optimizer_.set_option_value("verbose", get_bool_option("v"));
}

void
SequenceOptimizerApp::run() {

    // load motif graph from file
    auto lines =base::get_lines_from_file(get_string_option("mg"));
    auto mg = std::make_shared<motif_data_structure::MotifGraph>(lines[0],
                                                                 motif_data_structure::MotifGraphStringType::MG);

    // parse connection info
    _get_end_connections(mg);

    // get sequence optimizer scorer
    auto scorer = _setup_optimizer_scorer();

    for(auto & n : *mg) {
        if(n->data()->name().substr(0, 5) == "HELIX") {
            n->data()->mtype(util::MotifType::HELIX);
        }
    }

    mg->replace_ideal_helices();
    auto test_sols = optimizer_.get_optimized_sequences(mg, scorer);
    std::cout << test_sols[0]->dist_score << std::endl;

    auto dss = mg->designable_secondary_structure();
    mg->replace_helical_sequence(test_sols[0]->sequence);
    mg->write_pdbs("new");

    exit(0);

    std::ofstream out, sf_out;
    sf_out.open(get_string_option("score_file"));
    sf_out << "opt_num,opt_score,eterna_score,opt_sequence,opt_structure" << std::endl;
    
    out.open(get_string_option("out_file"));

    auto opt_num = 1;
    auto spl = base::split_str_by_delimiter(get_string_option("end_1"), ",");
    
    if(spl.size() != 2) {
        throw std::runtime_error(
            "incorrect format for end_1 must be graph_pos,end_pos, example 20,1 for the 20th node "
            "and end 1");
    }
    auto ni1 = std::stoi(spl[0]);
    auto ei1 = std::stoi(spl[1]);
    
    spl = base::split_str_by_delimiter(get_string_option("end_2"), ",");
    if(spl.size() != 2) {
        throw std::runtime_error(
            "incorrect format for end_2 must be graph_pos,end_pos, example 20,1 for the 20th node "
            "and end 1");
    }
    
    auto ni2 = std::stoi(spl[0]);
    auto ei2 = std::stoi(spl[1]);

    if(get_string_option("opt") == "Internal") {
        scorer = std::make_shared<sequence_optimization::InternalTargetScorer>(ni1, ei1, ni2, ei2, false);
    }
    
    auto mg_copy = std::make_shared<motif_data_structure::MotifGraph>(*mg);
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

void
SequenceOptimizerApp::_get_end_connections(
        motif_data_structure::MotifGraphOP mg) {
    connections_ = std::vector<ConnectionTemplate>();
    if(get_string_option("end_1") == "" && get_string_option("end_2") == "" && get_string_option("connections") == "") {
        throw SequenceOptimizerAppException("please supply either end_1 / end_2 or connections");
    }
    else if(get_string_option("end_1") != "" && get_string_option("end_2") != "") {
        connections_.push_back(_parse_end_commandline_args());
    }
    else if(get_string_option("connections") != "") {
        auto connection_strs = base::split_str_by_delimiter(get_string_option("connections"), ";");
        for(auto const & connection_str : connection_strs) {
            auto spl = base::split_str_by_delimiter(connection_str, " ");
            auto ni1 = std::stoi(spl[0]);
            auto ei1 = mg->get_node(ni1)->data()->get_end_index(spl[1]);

            auto ni2 = std::stoi(spl[2]);
            auto ei2 = mg->get_node(ni2)->data()->get_end_index(spl[3]);

            auto start = NodeIndexandEdge{ni1, ei1};
            auto end   = NodeIndexandEdge{ni2, ei2};
            connections_.push_back(ConnectionTemplate{start, end, "Internal"});
            break;
        }

        if(connections_.size() == 0) {
            throw SequenceOptimizerAppException("connections were supplied but were not processed correctly");
        }

    }
    else {
        throw SequenceOptimizerAppException("please supply either end_1 / end_2 or connections");
    }
}

ConnectionTemplate
SequenceOptimizerApp::_parse_end_commandline_args() {
    auto spl = base::split_str_by_delimiter(get_string_option("end_1"), ",");

    if(spl.size() != 2) {
        throw std::runtime_error(
                "incorrect format for end_1 must be graph_pos,end_pos, example 20,1 for the 20th node "
                        "and end 1");
    }
    auto ni1 = std::stoi(spl[0]);
    auto ei1 = std::stoi(spl[1]);

    spl = base::split_str_by_delimiter(get_string_option("end_2"), ",");
    if(spl.size() != 2) {
        throw std::runtime_error(
                "incorrect format for end_2 must be graph_pos,end_pos, example 20,1 for the 20th node "
                        "and end 1");
    }

    auto ni2 = std::stoi(spl[0]);
    auto ei2 = std::stoi(spl[1]);

    auto start = NodeIndexandEdge{ni1, ei1};
    auto end   = NodeIndexandEdge{ni2, ei2};
    return ConnectionTemplate{start, end, "Internal"};
}

sequence_optimization::SequenceOptimizerScorerOP
SequenceOptimizerApp::_setup_optimizer_scorer() {
    if(connections_.size() == 1) {
        auto c = connections_[0];
        return std::make_shared<sequence_optimization::InternalTargetScorer>(
                c.start.ni, c.start.ei, c.end.ni, c.end.ei, false);
    }
    else {
        auto sub_scorers = std::vector<sequence_optimization::SequenceOptimizerScorerOP>();
        for(auto const & c : connections_) {
            auto scorer = std::make_shared<sequence_optimization::InternalTargetScorer>(
                    c.start.ni, c.start.ei, c.end.ni, c.end.ei, false);
            sub_scorers.push_back(scorer);
        }
        return std::make_shared<sequence_optimization::MultiTargetScorer>(sub_scorers);
    }
}



int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);
    
    //load tectos
    auto tecto_dir = String(base::base_dir()+"/rnamake/lib/RNAMake/apps/simulate_tectos");
    resources::Manager::instance().add_motif(tecto_dir+"/resources/GAAA_tetraloop");
    resources::Manager::instance().add_motif(tecto_dir+"/resources/GGAA_tetraloop");

    
    auto app = SequenceOptimizerApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;

}
