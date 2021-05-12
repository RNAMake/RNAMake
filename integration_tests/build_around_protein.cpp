//
// Created by Hassan Abdelsamad on 4/28/21.
//

#include <CLI/CLI.hpp>

#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"
#include "design_rna_scaffold/design_rna_scaffold.cc"
#include <data_structure/graph_base.h>
#include "tools.h"


#include "common.hpp"


TEST_CASE ("build_around_protein") {
    auto base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    std::set_terminate(base::print_backtrace);
    auto app = DesignRNAScaffold();
    app.setup_options();
//    base::init_logging_with_file(app.log_level());
    auto args = "--pdb pp7_w_helix.pdb --start_bp L1-L25 --end_bp A116-A205 --log_level debug --sequences_per_design 1 --designs 5 --search_cutoff 15.0 --dump_pdbs --skip_sequence_optimization --search_max_size 75";
    app.app_.parse(args);
    app.run();

    SUBCASE ("Compare logs") {
        CHECK(check_logs("fasd", "dasd"));
    }
}