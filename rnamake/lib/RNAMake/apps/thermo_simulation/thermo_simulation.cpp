//
//  thermo_simulation.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "thermo_simulation.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "thermo_fluctuation/thermo_fluc_simulation_devel.h"

void
ThermoSimulationApp::setup_options() {
    add_option("mt", String(""), OptionType::STRING, true);
    add_option("n1", 0, OptionType::INT, true);
    add_option("e1", 0, OptionType::INT, true);
    add_option("n2", 0, OptionType::INT, true);
    add_option("e2", 0, OptionType::INT, true);
    add_option("extra_me", String(""), OptionType::STRING, false);
    add_option("start_pdbs", false, OptionType::BOOL, false);
}


void
ThermoSimulationApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    //cl_parser_.assign_options(cl_options _, optimizer_.options(), "optimizer");
}

void
ThermoSimulationApp::run() {
    if(get_string_option("extra_me") != "") {
        std::cout << "THERMO_SIMULATION: registered extra motif ensembles from file: ";
        std::cout << get_string_option("extra_me") << std::endl;
        RM::instance().register_extra_motif_ensembles(get_string_option("extra_me"));
    }

    auto lines =base::get_lines_from_file(get_string_option("mt"));
    auto mt = std::make_shared<MotifTree>(lines[0], MotifTreeStringType::MT_STR);
    auto mset = std::make_shared<MotifStateEnsembleTree>(mt);

    if(get_bool_option("start_pdbs")) {
        auto mt = mset->to_mst()->to_motif_tree();
        mt->write_pdbs();
        std::cout << "THERMO_SIMULATION: outputing each motif as nodes.*.pdb" << std::endl;

    }

    auto tfs = ThermoFlucSimulationDevel();
    tfs.set_option_value("steps", 1000000);
    tfs.setup(mset, get_int_option("n1"), get_int_option("n2"), get_int_option("e1"), get_int_option("e2"));
    auto score = tfs.run();
    std::cout << score << std::endl;

    
    
}

int main(int argc, const char * argv[]) {
    
    //load tectos
    auto tecto_dir = String(base_dir()+"/rnamake/lib/RNAMake/apps/simulate_tectos");
    RM::instance().add_motif(tecto_dir+"/resources/GAAA_tetraloop");
    RM::instance().add_motif(tecto_dir+"/resources/GGAA_tetraloop");
    
    auto app = ThermoSimulationApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;

}