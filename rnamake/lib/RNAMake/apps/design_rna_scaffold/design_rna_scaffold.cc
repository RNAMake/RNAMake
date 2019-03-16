//
// Created by Joseph Yesselman on 3/9/19.
//

#include "base/backtrace.hpp"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"

DesignRNAScaffold::DesignRNAScaffold():
        base::Application(),
        rm_(resources::Manager::instance()),
        parameters_(DesignRNAScaffoldParameters()) {

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
DesignRNAScaffold::setup_options() {
    // core inputs
    // from pdb
    add_option("pdb", String(""), base::OptionType::STRING, false);
    add_option("start_bp", String(""), base::OptionType::STRING, false);
    add_option("end_bp", String(""), base::OptionType::STRING, false);
    // from motif graph (not used much)
    add_option("mg", String(""), base::OptionType::STRING, false);

    // common options
    add_option("designs", 1, base::OptionType::INT, false);
    add_option("dump_pdbs", false, base::OptionType::BOOL, false);
    add_option("dump_scaffold_pdbs", false, base::OptionType::BOOL, false);
    add_option("helix_type", "ideal_helices", base::OptionType::STRING, false);
    add_option("search_type", "astar", base::OptionType::STRING, false);
    add_option("search_cutoff", 5.0f, base::OptionType::FLOAT, false);
    add_option("skip_sequence_optimization", false, base::OptionType::BOOL, false);
    add_option("sequence_opt_cutoff", 5.0f, base::OptionType::FLOAT, false);

    // less common options
    add_option("no_basepair_checks", false, base::OptionType::BOOL, false);

}

void
DesignRNAScaffold::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);

    parameters_.pdb    = get_string_option("pdb");    parameters_.start_bp = get_string_option("start_bp");
    parameters_.end_bp = get_string_option("end_bp"); parameters_.mg       = get_string_option("mg");
    parameters_.skip_sequence_optimization = get_bool_option("skip_sequence_optimization");
    parameters_.no_basepair_checks         = get_bool_option("no_basepair_checks");
}

void
DesignRNAScaffold::run() {
    if     (parameters_.pdb != "") { _setup_from_pdb(); }
    else if(parameters_.mg  != "") {}
    else                           {}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
DesignRNAScaffold::_setup_from_pdb() {
    auto struc = rm_.get_structure(parameters_.pdb, "scaffold");
    LOGI << "loaded pdb from file: " << parameters_.pdb;

    if (parameters_.start_bp == "" || parameters_.end_bp == "") {
        LOG_ERROR << "must supply the name of the start_bp and end_bp using option -start_bp and "
                    "-end_bp respectively when using -pdb option";
    }

    check_bp(parameters_.start_bp, struc, "start");
    check_bp(parameters_.end_bp,  struc, "end");

    mg_ = std::make_shared<motif_data_structure::MotifGraph>();

}

void
DesignRNAScaffold::check_bp(
        String const & name,
        structure::RNAStructureOP struc,
        String const & type) {


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(
        int argc,
        const char ** argv) {

    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    //start logging
    base::init_logging();

    //load extra motifs being used
    String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    auto app = DesignRNAScaffold();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}