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
#include "util/csv.h"

#include "common.hpp"


TEST_CASE ("build around protein") {
    auto base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    std::set_terminate(base::print_backtrace);
    auto app = DesignRNAScaffold();
    app.setup_options();
    base::init_logging(app.log_level());
    auto args = "--pdb min_tetraloop_receptor.pdb --start_bp A220-A253 --end_bp A144-A159 --log_level debug --sequences_per_design 1 --designs 3 --search_type mc --motif_path flex_helices,twoway,flex_helices,p4_p6_top,flex_helices,twoway,flex_helices --dump_pdbs --extra_pdbs p4_p6_top.pdb --search_cutoff 5 --max_helix_length 15";
    app.app_.parse(args);
    app.run();

    SUBCASE ("Comparing default") {
        auto *fp1 = fopen("orig_default.scores", "r");
        auto *fp2 = fopen("default.scores", "r");
        CHECK(compareFile(fp1, fp2));
    }

}