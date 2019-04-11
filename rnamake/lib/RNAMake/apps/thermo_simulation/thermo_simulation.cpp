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

#include <thermo_fluctuation/graph/simulation.h>

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

    add_option("log_level", "info", base::OptionType::STRING);

}


void
ThermoSimulationApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    base::Application::parse_command_line(argc, argv);

    parameters_.mg_file   = get_string_option("mg_file");
    parameters_.log_level = get_string_option("log_level");

}

void
ThermoSimulationApp::run() {
    // setup logging
    auto log_level = base::log_level_from_str(parameters_.log_level);
    base::init_logging(log_level);
    LOG_INFO << "log level set to: " << parameters_.log_level;

    auto lines = base::get_lines_from_file(parameters_.mg_file);
    LOG_INFO << lines.size()-1 << " designs loaded from " << parameters_.mg_file;

    auto mg = motif_data_structure::MotifGraphOP(nullptr);
    auto scorer = std::make_shared<thermo_fluctuation::graph::FrameScorer>();
    auto sim = std::make_shared<thermo_fluctuation::graph::Simulation>(scorer);
    for(auto const & l : lines) {
        // make sure its not the dummy line at the end
        if(l.size() < 10) { break; }

        mg = std::make_shared<motif_data_structure::MotifGraph>(l, motif_data_structure::MotifGraphStringType::MG);
        auto r = _get_mseg(mg);
        //std::cout << mseg->size() << " " << mg->size() << std::endl;
        auto c = _guess_connection(mg);
        std::cout << r->index_hash[c.start.node_index] << " " << r->index_hash[c.end.node_index] << std::endl;

        auto start = data_structure::NodeIndexandEdge{r->index_hash[c.start.node_index], c.start.edge_index};
        auto end = data_structure::NodeIndexandEdge{r->index_hash[c.end.node_index], c.end.edge_index};

        int avg = 0;
        for(int j = 0; j < 5; j++) {
            sim->setup(*r->mseg, start, end);
            int count = 0;
            for (int i = 0; i < 1000000; i++) {
                count += sim->next();
            }
            avg += count;
        }
        avg /= 5;
        std::cout << avg << std::endl;


    }

    exit(0);

    if(get_string_option("extra_me") != "") {
        std::cout << "THERMO_SIMULATION: registered extra motif ensembles from file: ";
        std::cout << get_string_option("extra_me") << std::endl;
        resources::Manager::instance().register_extra_motif_ensembles(get_string_option("extra_me"));
    }

    /*auto lines =base::get_lines_from_file(get_string_option("mt"));
    auto mt = std::make_shared<motif_data_structure::MotifTree>(lines[0], motif_data_structure::MotifTreeStringType::MT_STR);
    auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>(mt);

    if(get_bool_option("start_pdbs")) {
        auto mt = mset->to_mst()->to_motif_tree();
        mt->write_pdbs();
        std::cout << "THERMO_SIMULATION: outputing each motif as nodes.*.pdb" << std::endl;

    }

    auto tfs = thermo_fluctuation::ThermoFlucSimulationDevel();
    tfs.set_option_value("steps", 1000000);
    tfs.setup(mset, get_int_option("n1"), get_int_option("n2"), get_int_option("e1"), get_int_option("e2"));
    auto score = tfs.run();
    std::cout << score << std::endl;*/

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
            int pi = index_hash[n->parent_index()];
            j = mseg->add_ensemble(*mse, data_structure::NodeIndexandEdge{pi, n->parent_end_index()});
        }

        index_hash[n->index()] = j;

    }
    return std::make_shared<EnsembleConversionResults>(mseg, index_hash);

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