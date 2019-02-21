//
// Created by Joseph Yesselman on 2/16/19.
//

#ifndef TEST_DOCK_MOTIF_H
#define TEST_DOCK_MOTIF_H

#include <stdio.h>
#include "base/application.hpp"

class DockMotifAppException : public std::runtime_error {
public:
    DockMotifAppException(
            String const & message):
            std::runtime_error(message)
    {}

};

class DockMotifApp : public Application {
public:
    DockMotifApp();

    ~DockMotifApp() {}

public: // application interface functions

    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:
    ResidueOP
    _parse_ligand_for_center_coords();

    Point
    _calc_motif_center(
            MotifOP);

};


#endif //TEST_DOCK_MOTIF_H
