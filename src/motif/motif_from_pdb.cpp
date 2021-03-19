#include "motif/motif_from_pdb.h"
#include "motif/motif_factory.h"
#include "motif"
#include <vector>
#include <filesystem>

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;

namespace motif {

    MotifOPs motifs_;

void
compute_motif(String const & path, vector<int> list) {
    // Add a for loop
    for (const auto& file_path : recursive_directory_iterator(path)) {
        MotifFactory factory = MotifFactory();
        auto motif = factory.motif_from_file(file_path);
        motif = factory.align_motif_to_common_frame(motif, 0);
        motifs_.push_back(motif);
    }

    // generate ensemble with with a list of energies
    // use to_str for list
    auto motifEnsemble = MotifEnsemble(motifs_, list);
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