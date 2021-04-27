#include <CLI/CLI.hpp>

#include "motif_from_pdb.h"
#include "motif/motif_factory.h"
#include "motif/motif_ensemble.h"
#include "motif/motif.h"
#include "base/backtrace.h"
#include "base/log.h"
#include <CLI/CLI.hpp>

#include <vector>
#include <filesystem>

using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;


ComputeEnsemble::ComputeEnsemble() :
        base::Application(),
        rm_(resources::Manager::instance()),
        app_("ComputeEnsemble") {}

String
valid_pdb(String &path) {
    auto ending = path.substr(path.size() - 4);
    return ending == ".pdb" ? String{""} : String{"the file specified by --pdb must end in .pdb"};
}

void
ComputeEnsemble::compute_motif(String const &path, std::vector<float> const &list) {
    // Add a for loop
    for (const auto &file_path_entry : recursive_directory_iterator(path)) {
        auto path_string = file_path_entry.path().string();
        auto factory = motif::MotifFactory();
        auto motif = factory.motif_from_file(path_string);
        motif = factory.align_motif_to_common_frame(motif, 0);
        motifs_.push_back(motif);
    }

    // generate ensemble with with a list of energies
    // use to_str for list
    auto motif_ensemble = motif::MotifEnsemble(motifs_, list);
    auto s = motif_ensemble.to_str();
    std::fstream file;
    file.open("test.txt");
    file << s;
    file.close();
    std::cout << s << std::endl;
}

void
ComputeEnsemble::run() {
    std::vector<float> vec;
    vec.push_back(10.0);
    vec.push_back(15.0);
    ComputeEnsemble::compute_motif(parameters_.pdb, vec);
}

void
ComputeEnsemble::parse_command_line (
        int argc,
        const char **argv) {
}

void
ComputeEnsemble::setup_options () {
    app_.add_option("--pdb", parameters_.pdb, "path to a PDB file with input RNA structure");
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
