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

    std::vector<motif::MotifStateOPs>
    _get_libraries(
            String const &);


private:
    util::StericLookup lookup_;
};


#endif //TEST_APT_NEW_INTERFACE_H
