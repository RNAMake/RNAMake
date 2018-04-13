//
// Created by Joseph Yesselman on 4/12/18.
//

#ifndef TEST_APT_STABLIZATION_H
#define TEST_APT_STABLIZATION_H

#include <stdio.h>
#include "base/application.hpp"
#include "util/steric_lookup.hpp"

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

public:

};

#endif //TEST_APT_STABLIZATION_H
