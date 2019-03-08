//
// Created by Joseph Yesselman on 3/7/19.
//

#ifndef TEST_SEQUENCE_OPTIMIZATION_BENCHMARKS_H
#define TEST_SEQUENCE_OPTIMIZATION_BENCHMARKS_H

#include "base/application.hpp"
#include "util/steric_lookup.hpp"
#include "motif_data_structures/motif_state_graph.hpp"
#include "motif_state_search/motif_state_monte_carlo.h"

class SequenceOptProblem {
public:
    virtual
    MotifStateGraphOP
    get_motif_state_graph();



};





class SequenceOptimizationBenchmarks : public Application {
public:
    SequenceOptimizationBenchmarks();

    ~SequenceOptimizationBenchmarks() {}

public:
    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:
    MotifStateGraphOP
    _get_starting_graph();

    void
    _setup_sterics(
            MotifStateGraphOP,
            StericLookup &);

    std::vector<MotifStateOPs>
    _get_libraries(
            String const &);


};


#endif //TEST_SEQUENCE_OPTIMIZATION_BENCHMARKS_H
