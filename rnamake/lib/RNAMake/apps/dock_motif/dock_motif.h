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

private: // search functions
    void
    _score(
            MotifStateOP);

    void
    _search();

    MotifStateOP
    _get_starting_state(
            StericLookup &,
            Point const &);

private:
    ResidueOP
    _parse_ligand_for_center_coords();

    Point
    _calc_motif_center(
            MotifOP);

    Matrix
    _rotation_about_x_axis(
            float);

    Matrix
    _rotation_about_y_axis(
            float);

    Matrix
    _rotation_about_z_axis(
            float);

    void
    _precompute_rotations(
            MotifStateOP);

private:
    MotifStateOPs rotations_;
    StericLookup lookup_;
    Point center_;

};

Point
get_random_point(
        RandomNumberGenerator &,
        int);

#endif //TEST_DOCK_MOTIF_H
