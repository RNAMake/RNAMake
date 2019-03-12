//
// Created by Joseph Yesselman on 3/9/19.
//

#ifndef TEST_DESIGN_RNA_SCAFFOLD_H
#define TEST_DESIGN_RNA_SCAFFOLD_H

#include <stdio.h>
#include "base/application.hpp"
#include "util/steric_lookup.hpp"
#include "motif_state_search/motif_state_search.h"
#include "motif_data_structures/motif_graph.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"

class DesignRNAScaffoldAppException : public std::runtime_error {
public:
    DesignRNAScaffoldAppException(
            String const & message):
            std::runtime_error(message)
    {}
};


struct EndStateInfo {
    String name;
    int n_pos;
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

    std::shared_ptr<MSS_Path>
    _setup_path();

    std::vector<MotifStateOPs>
    _get_libraries();

private:
    MotifStateSearch search_;
    MotifGraphOP mg_;
    EndStateInfo start_, end_;
    util::StericLookup lookup_;
    SequenceOptimizer3D optimizer_;


};



#endif //TEST_DESIGN_RNA_SCAFFOLD_H
