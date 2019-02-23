//
// Created by Joseph Yesselman on 2/23/19.
//

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_graph.hpp"
#include "motif_state_search/motif_state_monte_carlo.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"
#include "apt_new_interface/apt_new_interface.h"


AptNewInterface::AptNewInterface() {}

void
AptNewInterface::setup_options() {
    add_option("scaffold", String(""), OptionType::STRING, true);
    add_option("docked_motif", String(""), OptionType::STRING, true);

    // general options
    add_option("out_file", "default.out", OptionType::STRING, false);
    add_option("score_file", "default.scores", OptionType::STRING, false);
    add_option("designs", 1, OptionType::INT, false);

}

void
AptNewInterface::parse_command_line(
        int argc,
        const char **argv) {

    Application::parse_command_line(argc, argv);
}

void
AptNewInterface::run() {

    auto start_path = String("flex_helices,twoway,flex_helices");



}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(
        int argc,
        const char ** argv) {

    //must add this for all apps!
    std::set_terminate(print_backtrace);

    //load extra motifs being used
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/apt_new_interface/resources/";
    RM::instance().add_motif(base_path+"pRNA_3WJ");

    auto app = AptNewInterface();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}