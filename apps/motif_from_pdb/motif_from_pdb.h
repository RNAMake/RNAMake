#include <stdio.h>
#include <regex>

#include "base/application.hpp"
#include <data_structure/graph_base.h>
#include "base/application.hpp"
#include "util/steric_lookup.hpp"
#include "motif_search/search.h"
#include "motif_search/solution_topology.h"
#include <motif_search/solution_filter.h>
#include "motif_data_structure/motif_graph.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"
#include <thermo_fluctuation/graph/simulation.h>

class ComputeEnsemble : public base::Application {

public: // application functions

    void
    setup_options() override;

    void
    parse_command_line(
            int,
            const char **) override;

    void
    run() override;

    void
    compute_motif(String const &path, std::vector<float> const &list);

public:

    ComputeEnsemble();

    CLI::App app_;

private:


    struct Parameters {
        String pdb = "";
    };

private:
    Parameters parameters_ = Parameters();
    resources::Manager &rm_;
    motif::MotifOPs motifs_;


};