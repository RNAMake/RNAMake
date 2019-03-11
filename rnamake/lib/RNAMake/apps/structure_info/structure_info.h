//
// Created by Joseph Yesselman on 11/21/18.
//

#ifndef TEST_STRUCTURE_INFO_H
#define TEST_STRUCTURE_INFO_H

#include <stdio.h>

//RNAMake Headers
#include "base/application.hpp"

class StructureInfoAppException : public std::runtime_error {
public:
    StructureInfoAppException(
            String const & message):
            std::runtime_error(message) {}
};

class StructureInfoApp : public base::Application {
public:
    StructureInfoApp();

public: // application setups functions
    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

public:

    void
    run();

};


#endif //TEST_STRUCTURE_INFO_H
