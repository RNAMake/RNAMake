//
// Created by Hassan Abdelsamad on 4/18/21.
//

#include <CLI/CLI.hpp>

#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"
#include "design_rna_scaffold/design_rna_scaffold.cc"
#include <data_structure/graph_base.h>
#include "tools.h"

#include "common.hpp"


TEST_CASE ("Test TTR") {
    auto base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    std::set_terminate(base::print_backtrace);
    auto app = DesignRNAScaffold();

    app.setup_options();
    base::init_logging(app.log_level());
    auto args = "--pdb min_tetraloop_receptor.pdb --start_bp A220-A253 --end_bp A144-A159 --log_level debug --sequences_per_design 1 --designs 10 --search_cutoff 5.0";
    app.app_.parse(args);
    app.run();

    SUBCASE ("Compare default.scores") {
        auto *fp1 = fopen("orig_default.scores", "r");
        auto *fp2 = fopen("default.scores", "r");
        CHECK(compareFile(fp1, fp2));
    }

    SUBCASE ("Compare default.out") {
        auto *fp1 = fopen("orig_default.out", "r");
        auto *fp2 = fopen("default.out", "r");
        CHECK(compareFile(fp1, fp2));
    }
}

