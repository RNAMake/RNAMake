//
// Created by Hassan Abdelsamad on 4/28/21.
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

TEST_CASE ("build_around_protein") {
    String base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    std::set_terminate(base::print_backtrace);
    auto app = DesignRNAScaffold();
    app.setup_options();
    base::init_logging(app.log_level());
    auto args = "--pdb pp7_w_helix.pdb --start_bp L1-L25 --end_bp A116-A205 --log_level debug --sequences_per_design 1 --designs 5 --search_cutoff 15.0 --dump_pdbs --skip_sequence_optimization --search_max_size 75";
    app.app_.parse(args);
    app.run();

    SUBCASE ("Compare default.scores") {
        auto *fp1 = fopen("orig_default.scores", "r");
        auto *fp2 = fopen("default.scores", "r");
        CHECK(compareFile(fp1, fp2));

    }
}