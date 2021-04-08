#include <CLI/CLI.hpp>

#include "motif/motif_from_pdb.h"
#include "motif/motif_factory.h"
#include "motif/motif_ensemble.h"
#include "motif/motif.h"
#include "base/backtrace.h"
#include "base/log.h"
#include <CLI/CLI.hpp>

#include <CLI/CLI.hpp>

#include <vector>
#include <filesystem>

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;

motif::MotifOPs motifs_;

ComputeEnsemble::ComputeEnsemble() :
        base::Application(),
        rm_(resources::Manager::instance()),
        app_("ComputeEnsemble") {}

void
compute_motif(String const &path, std::vector<float> list) {
    // Add a for loop
    for (const auto &file_path_entry : recursive_directory_iterator(path)) {
        String path_string = file_path_entry.path().string();
        motif::MotifFactory factory = motif::MotifFactory();
        auto motif = factory.motif_from_file(path_string);
        motif = factory.align_motif_to_common_frame(motif, 0);
        motifs_.push_back(motif);
    }

    // generate ensemble with with a list of energies
    // use to_str for list
    auto motifEnsemble = motif::MotifEnsemble(motifs_, list);
}

int
main(int argc, const char **argv) {
    std::set_terminate(base::print_backtrace);
    auto app = ComputeEnsemble();
    app.setup_options();
    CLI11_PARSE(app.app_, argc, argv);
    app.run();

    return 0;

}