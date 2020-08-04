//
// Created by Joseph Yesselman on 2/23/19.
//

#ifndef TEST_APT_NEW_INTERFACE_H
#define TEST_APT_NEW_INTERFACE_H

#include <stdio.h>
#include "base/application.hpp"
#include "util/steric_lookup.hpp"
#include "motif_data_structure/motif_state_graph.hpp"

class AptNewInterfaceException : public std::runtime_error {
public:
    AptNewInterfaceException(
            String const & message):
            std::runtime_error(message)
    {}
};

struct NodeIndexandEdge {
    int ni; //node index
    int ei; //end index
};

struct ConnectionInfo {
    NodeIndexandEdge start;
    NodeIndexandEdge end;
};

struct AptNewInterfaceProblem {
    inline
    AptNewInterfaceProblem(
            motif_data_structure::MotifStateGraphOP n_msg,
            ConnectionInfo const & n_path_1,
            ConnectionInfo const & n_path_2):
            msg(n_msg),
            path_1(n_path_1),
            path_2(n_path_2) {}

    motif_data_structure::MotifStateGraphOP msg;
    ConnectionInfo path_1, path_2;

};

typedef std::shared_ptr<AptNewInterfaceProblem> AptNewInterfaceProblemOP;

class AptNewInterface : public base::Application {
public:
    AptNewInterface();

    ~AptNewInterface() {}

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
    void
    _setup_sterics(
            motif_data_structure::MotifStateGraphOP);

private:

    void
    _setup_new_motifs();

    AptNewInterfaceProblemOP
    _get_design_problem();

    motif_data_structure::MotifStateGraphOP
    _setup_graph();

    std::vector<motif::MotifStateOPs>
    _get_libraries(
            int);

    motif_search::MotifStateMonteCarloOP
    _setup_search(
            ConnectionInfo const &,
            util::StericLookupNew const &,
            motif_data_structure::MotifStateGraphOP,
            std::vector<motif::MotifStateOPs> const &,
            float);

    int
    _find_min_motifs_per_path(
            ConnectionInfo const &,
            util::StericLookupNew const &,
            motif_data_structure::MotifStateGraphOP,
            float,
            int);

    int
    _get_residue_count(
            motif_data_structure::MotifGraphOP);

    String
    _get_motif_names(
            motif_data_structure::MotifGraphOP);


private:
    util::StericLookupNew lookup_;
    util::RandomNumberGenerator rng_;
    resources::Manager & rm_;
};


#endif //TEST_APT_NEW_INTERFACE_H
