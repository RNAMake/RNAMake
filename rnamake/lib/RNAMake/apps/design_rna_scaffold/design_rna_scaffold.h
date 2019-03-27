//
// Created by Joseph Yesselman on 3/9/19.
//

#ifndef TEST_DESIGN_RNA_SCAFFOLD_H
#define TEST_DESIGN_RNA_SCAFFOLD_H

#include <stdio.h>
#include "base/application.hpp"
#include "util/steric_lookup.hpp"
#include "motif_search/path_finding/search.h"
#include "motif_data_structure/motif_graph.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"

struct EndStateInfo {
    String name;
    int n_pos;
};

struct Parameters {
    String pdb, start_bp, end_bp, mg;
    bool skip_sequence_optimization, no_basepair_checks;
};


class DesignRNAScaffold : public base::Application {
public:

    DesignRNAScaffold();

public: // application functions

    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:

    void
    _setup_sterics();

    void
    _setup_from_pdb();

    void
    _setup_from_mg();

    std::vector<motif::MotifStateOPs>
    _get_libraries();

private:
    void
    check_bp(
            String const &,
            structure::RNAStructureOP,
            String const &);

private:
    motif_search::SearchOP search_;
    motif_data_structure::MotifGraphOP mg_;
    EndStateInfo start_, end_;
    util::StericLookup lookup_;
    sequence_optimization::SequenceOptimizer3D optimizer_;
    resources::Manager & rm_;
    Parameters parameters_;


};



#endif //TEST_DESIGN_RNA_SCAFFOLD_H
