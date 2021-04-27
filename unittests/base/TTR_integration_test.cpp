//
// Created by Hassan Abdelsamad on 4/18/21.
//

#include <CLI/CLI.hpp>

#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"
#include "design_rna_scaffold/design_rna_scaffold.cc"
#include <data_structure/graph_base.h>

#include "../common.hpp"

TEST_CASE ("Test application") {
            SUBCASE("Test") {
        std::cout << "Hello" << "\n";
        std::set_terminate(base::print_backtrace);
        auto app = DesignRNAScaffold();

        app.setup_options();
//        base::init_logging(app.log_level());

        auto args = "--pdb inputs/min_tetraloop_receptor.pdb --start_bp A220-A253 --end_bp A144-A159 --log_level debug --sequences_per_design 1 --designs 10 --search_cutoff 5.0";
        app.app_.parse(args);

        app.run();
    }
}