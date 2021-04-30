//
// Created by Hassan Abdelsamad on 4/30/21.
//

#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"
#include "design_rna_scaffold/design_rna_scaffold.cc"
#include <data_structure/graph_base.h>


#include "common.hpp"

int compareFile(FILE *fPtr1, FILE *fPtr2) {
    char ch1, ch2;
    auto line = 1;
    auto col = 0;

    do {
        // Input character from both files
        ch1 = fgetc(fPtr1);
        ch2 = fgetc(fPtr2);

        // Increment line
        if (ch1 == '\n') {
            line += 1;
            col = 0;
        }

        // If characters are not same then return false
        if (ch1 != ch2)
            return -1;
        col += 1;

    } while (ch1 != EOF && ch2 != EOF);

    /* If both files have reached end */
    if (ch1 == EOF && ch2 == EOF)
        return 1;
    else
        return -1;
}


TEST_CASE ("Test TTR_MC") {
    auto base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    std::set_terminate(base::print_backtrace);
    auto app = DesignRNAScaffold();

    app.setup_options();
    base::init_logging(app.log_level());
    auto args = "--pdb inputs/min_tetraloop_receptor.pdb --start_bp A220-A253 --end_bp A144-A159 --log_level debug --sequences_per_design 1 --designs 3 --search_type mc --motif_path flex_helices,twoway,flex_helices,twoway,flex_helices  --dump_pdbs --thermo_fluc ";
    app.app_.parse(args);
    app.run();

        SUBCASE ("Compare default.scores") {
            auto *fp1 = fopen("orig_default.scores", "r");
            auto *fp2 = fopen("default.scores", "r");
            CHECK(compareFile(fp1, fp2));

        SUBCASE ("Compare default.out") {
            auto *fp1 = fopen("orig_default.out", "r");
            auto *fp2 = fopen("default.out", "r");
            CHECK(compareFile(fp1, fp2));
        }

        SUBCASE ("Compare design-0") {

        }
        SUBCASE ("Compare design-1") {

        }
        SUBCASE ("Compare design-2") {

        }
    }

