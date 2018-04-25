//
// Created by Joseph Yesselman on 4/12/18.
//

#ifndef TEST_APT_STABLIZATION_H
#define TEST_APT_STABLIZATION_H

#include <stdio.h>
#include "base/application.hpp"
#include "util/steric_lookup.hpp"

class APTStablizationException : public std::runtime_error {
public:
    APTStablizationException(
            String const & message):
            std::runtime_error(message)
    {}
};


class APTStablization : public Application {
public:
    APTStablization();

    ~APTStablization() {}

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
            MotifStateGraphOP);

private:

    std::vector<MotifStateOPs>
    _get_libraries(
            String const &);

private:
    StericLookup lookup_;

};

#endif //TEST_APT_STABLIZATION_H
