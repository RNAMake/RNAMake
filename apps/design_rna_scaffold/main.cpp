//
// Created by Hassan Abdelsamad on 4/18/21.
//

#include <CLI/CLI.hpp>

#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"


int
main (int argc, const char **argv) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    String base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");
    auto app = DesignRNAScaffold();
    app.setup_options();
    CLI11_PARSE(app.app_, argc, argv);
    //start logging
    base::init_logging(base::LogLevel::DEBUG);

    // hacky way of doing it but wtv, the app is guaranteed to have an input CJ 09/20
    //if(app.app_["--mg"]->empty() && app.app_["--pdb"]->empty()) {
    //    LOGF<<"you must input a PDB or motif graph file via --pdb or --mg";
    //    exit(1);
    //}

    //load extra motifs being used
    app.run();
    return 0;
}