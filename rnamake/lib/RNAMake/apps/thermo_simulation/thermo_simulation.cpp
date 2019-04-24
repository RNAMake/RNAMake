//
//  thermo_simulation.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "thermo_simulation.hpp"

#include <base/backtrace.hpp>
#include <base/log.h>

#include "base/file_io.h"
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_state_ensemble_tree.h"
#include "thermo_fluctuation/thermo_fluc_simulation_devel.h"

void
ThermoSimulationApp::setup_options() {
    add_option("mg_file", String(""), base::OptionType::STRING, true);
    /*add_option("mt", String(""), base::OptionType::STRING, true);
    add_option("n1", 0, base::OptionType::INT, true);
    add_option("e1", 0, base::OptionType::INT, true);
    add_option("n2", 0, base::OptionType::INT, true);
    add_option("e2", 0, base::OptionType::INT, true);
    add_option("extra_me", String(""), base::OptionType::STRING, false);
    add_option("start_pdbs", false, base::OptionType::BOOL, false);*/

    add_option("steps", 100000, base::OptionType::INT);
    add_option("n", 1, base::OptionType::INT);
    add_option("log_level", "info", base::OptionType::STRING);
    add_option("extra_sequences", "", base::OptionType::STRING);
    add_option("score_file", "thermo_sim.scores", base::OptionType::STRING);
}


void
ThermoSimulationApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    base::Application::parse_command_line(argc, argv);

    parameters_.mg_file    = get_string_option("mg_file");
    parameters_.log_level  = get_string_option("log_level");
    parameters_.steps      = get_int_option("steps");
    parameters_.n          = get_int_option("n");
    parameters_.extra_sequences = get_string_option("extra_sequences");
    parameters_.score_file = get_string_option("score_file");

}

void
ThermoSimulationApp::run() {
    // setup logging
    auto log_level = base::log_level_from_str(parameters_.log_level);
    base::init_logging(log_level);
    LOG_INFO << "log level set to: " << parameters_.log_level;

    if(parameters_.extra_sequences != "") { _parse_extra_sequences(parameters_.extra_sequences); }
    auto lines = base::get_lines_from_file(parameters_.mg_file);
    LOG_INFO << lines.size()-1 << " designs loaded from " << parameters_.mg_file;

    auto score_out = std::ofstream();
    score_out.open(parameters_.score_file);
    score_out << "design_num,hits,opt_score,opt_sequence,opt_structure" << std::endl;

    auto mg = motif_data_structure::MotifGraphOP(nullptr);
    auto scorer = std::make_shared<thermo_fluctuation::graph::FrameScorer>();
    auto sterics = std::make_shared<thermo_fluctuation::graph::sterics::NoSterics>();
    auto sim = std::make_shared<thermo_fluctuation::graph::Simulation>(scorer, sterics);
    scorer->setup(true);
    int design_count = 0;
    for(auto const & l : lines) {
        // make sure its not the dummy line at the end
        if(l.size() < 10) { break; }

        mg = std::make_shared<motif_data_structure::MotifGraph>(l, motif_data_structure::MotifGraphStringType::MG);
        auto r = _get_mseg(mg);
        auto c = _guess_connection(mg);


        auto start = data_structure::NodeIndexandEdge{r->index_hash[c.start.node_index], c.start.edge_index};
        auto end = data_structure::NodeIndexandEdge{r->index_hash[c.end.node_index], c.end.edge_index};

        int avg = _calc_num_hits(sim, *r->mseg, start, end, parameters_.n);
        LOG_INFO << "design num: " << design_count << " hits target " << avg;
        score_out << design_count << "," << avg << "," << _calc_initial_score(scorer, mg, c.start, c.end) << ",";
        score_out << mg->sequence() << "," << mg->dot_bracket() << std::endl;
        if(extra_sequences_.find(design_count) == extra_sequences_.end()) {
            design_count += 1;
            continue;
        }
        for(auto const & seq : extra_sequences_[design_count]) {
            mg->replace_helical_sequence(seq);
            int avg = _calc_num_hits(sim, *r->mseg, start, end, parameters_.n);
            LOG_INFO << "design num: " << design_count << " hits target " << avg << " with sequence: " << seq;
            score_out << design_count << "," << avg << "," << _calc_initial_score(scorer, mg, c.start, c.end) << ",";
            score_out << mg->sequence() << "," << mg->dot_bracket() << std::endl;
        }

        design_count += 1;
    }

}

ConnectionTemplate
ThermoSimulationApp::_guess_connection(
        motif_data_structure::MotifGraphOP mg) {

    auto connections = std::vector<ConnectionTemplate>();
    for(auto const & n : *mg) {
        for(auto const & c : n->connections()) {
            if(c == nullptr) { continue; }

            auto ei1 = c->end_index(n->index());
            auto partner = c->partner(n->index());
            auto ei2 = c->end_index(partner->index());

            auto end1 = n->data()->ends()[ei1];
            auto end2 = partner->data()->ends()[ei2];
            auto diff = end1->state()->diff(end2->state());

            if(diff < 3) { continue; }

            auto nie1 = data_structure::NodeIndexandEdge{n->index(), ei1};
            auto nie2 = data_structure::NodeIndexandEdge{partner->index(), ei2};

            if(n->index() < partner->index()) { std::swap(nie1, nie2); }

            auto ct = ConnectionTemplate{nie1, nie2};
            bool found = false;
            for(auto const & ct2 : connections) {
                if(ct == ct2) { found = true; }
            }
            if(!found) {
                connections.push_back(ct);
            }
        }
    }

    if(connections.size() > 1) {
        LOG_ERROR << "multiple connections possible, please specify at commandline"; exit(0);
    }

    if(connections.size() == 0) {
        LOG_ERROR << "no connections possible, please specify or try another MotifGraph"; exit(0);
    }

    return connections[0];
}


ThermoSimulationApp::EnsembleConversionResultsOP
ThermoSimulationApp::_get_mseg(
        motif_data_structure::MotifGraphOP mg) {

    auto index_hash = std::map<int, int>();
    auto mseg = std::make_shared<motif_data_structure::MotifStateEnsembleGraph>();
    int i = 0, j = 0;
    for(auto const & n : *mg) {
        // build / get motif state ensemble
        auto mse = motif::MotifStateEnsembleOP(nullptr);

        if(n->data()->mtype() == util::MotifType::HELIX) {
            // is not a basepair step
            if(n->data()->residues().size() > 4) {
                LOG_ERROR << "supplied a helix motif: "+n->data()->name() +" that is not a basepair step, this is no supported";
                exit(0);
            }

            try {
                mse = rm_.motif_state_ensemble(n->data()->end_ids()[0]);
            }
            catch(resources::ResourceManagerException const & e) {
                LOG_ERROR << "cannot find motif state ensemble for basepair with id: " + n->data()->end_ids()[0] <<
                "check to make sure its a Watson-Crick basepair step";
                exit(0);
            }
        }
        else {
            mse = std::make_shared<motif::MotifStateEnsemble>(n->data()->get_state());
        }

        if(i == 0) {
            j = mseg->add_ensemble(*mse);
        }
        else {
            int pi = index_hash[mg->parent_index(n->index())];
            int pie = mg->parent_end_index(n->index());
            j = mseg->add_ensemble(*mse, data_structure::NodeIndexandEdge{pi, pie});
        }

        index_hash[n->index()] = j;
        i++;

    }
    return std::make_shared<EnsembleConversionResults>(mseg, index_hash);

}

void
ThermoSimulationApp::_parse_extra_sequences(
        String const & file_path) {
    extra_sequences_ = std::map<int, Strings>();
    auto lines = base::get_lines_from_file(file_path);
    LOG_INFO << "loading extra sequences from: " << file_path << " contains " << lines.size()-2 << " sequences";
    auto spl = Strings();
    int i = -1;
    int design_num;
    for(auto const & l : lines) {
        i++;
        // make sure its not the dummy line at the end
        if(l.size() < 10) { break; }
        if(i == 0) { continue; }
        spl = base::split_str_by_delimiter(l, ",");
        design_num = std::stoi(spl[0]);
        if(extra_sequences_.find(design_num) == extra_sequences_.end()) {
            extra_sequences_[design_num] = Strings();
        }
        extra_sequences_[design_num].push_back(spl[1]);
    }
}

int
ThermoSimulationApp::_calc_num_hits(
        thermo_fluctuation::graph::SimulationOP sim,
        motif_data_structure::MotifStateEnsembleGraph const & mseg,
        data_structure::NodeIndexandEdge const & start,
        data_structure::NodeIndexandEdge const & end,
        int n) {
    int avg = 0;
    for(int j = 0; j < n; j++) {
        sim->setup(mseg, start, end);
        int count = 0;
        for (int i = 0; i < parameters_.steps; i++) {
            count += sim->next();
        }
        avg += count;
    }
    avg /= n;
    return avg;
}

float
ThermoSimulationApp::_calc_initial_score(
        thermo_fluctuation::graph::ScorerOP scorer,
        motif_data_structure::MotifGraphOP mg,
        data_structure::NodeIndexandEdge const & start,
        data_structure::NodeIndexandEdge const & end) {
    auto end_state_1 = mg->get_node(start.node_index)->data()->ends()[start.edge_index]->state();
    auto end_state_2 = mg->get_node(end.node_index)->data()->ends()[end.edge_index]->state();
    return scorer->score(*end_state_1, *end_state_2);
}



int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);
    
    //load tectos
    auto tecto_dir = String(base::base_dir()+"/rnamake/lib/RNAMake/apps/simulate_tectos");
    resources::Manager::instance().add_motif(tecto_dir+"/resources/GAAA_tetraloop");
    resources::Manager::instance().add_motif(tecto_dir+"/resources/GGAA_tetraloop");
    
    auto app = ThermoSimulationApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;

}