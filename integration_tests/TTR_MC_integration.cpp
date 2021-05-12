//
// Created by Hassan Abdelsamad on 4/30/21.
//

#include <CLI/CLI.hpp>

#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"
#include "design_rna_scaffold/design_rna_scaffold.cc"
#include <data_structure/graph_base.h>
#include "tools.h"


#include "common.hpp"


TEST_CASE ("Test TTR_MC") {
    auto base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    std::set_terminate(base::print_backtrace);
    auto app = DesignRNAScaffold();

    app.setup_options();
    base::init_logging_with_file(app.log_level());
    auto args = "--pdb min_tetraloop_receptor.pdb --start_bp A220-A253 --end_bp A144-A159 --log_level debug --sequences_per_design 1 --designs 3 --search_type mc --motif_path flex_helices,twoway,flex_helices,twoway,flex_helices  --dump_pdbs --thermo_fluc ";
    app.app_.parse(args);
    app.run();

        SUBCASE ("Compare default.scores") {

        }

        SUBCASE ("Compare default.out") {

        }

        SUBCASE ("Compare design-0") {

        }
        SUBCASE ("Compare design-1") {

        }
        SUBCASE ("Compare design-2") {

        }
    }
