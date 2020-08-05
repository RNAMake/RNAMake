//
// Created by Joseph Yesselman on 2/16/19.
//

#ifndef TEST_DOCK_MOTIF_H
#define TEST_DOCK_MOTIF_H

#include <stdio.h>
#include "base/application.hpp"

struct MotifStateandScore {
    motif::MotifStateOP ms;
    float score;
};

bool sort_by_score(
        MotifStateandScore const & a,
        MotifStateandScore const & b) {
    return a.score < b.score;
}

typedef std::vector<MotifStateandScore> MotifStateandScores;


class DockMotifAppException : public std::runtime_error {
public:
    DockMotifAppException(
            String const & message):
            std::runtime_error(message)
    {}

};

class DockMotifApp : public base::Application {
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
    float
    _score(
            motif::MotifStateOP);

    void
    _search();

    motif::MotifStateOP
    _get_starting_state(
            util::StericLookupNew &,
            math::Point const &);

private:
    structure::ResidueOP
    _parse_ligand_for_center_coords();

    math::Point
    _calc_motif_center(
            motif::MotifOP);

    math::Matrix
    _rotation_about_x_axis(
            float);

    math::Matrix
    _rotation_about_y_axis(
            float);

    math::Matrix
    _rotation_about_z_axis(
            float);

    void
    _precompute_rotations(
            motif::MotifStateOP);

private:
    motif::MotifStateOPs rotations_;
    util::StericLookupNew lookup_;
    math::Point center_;
    MotifStateandScores results_;
    motif::MotifStateOP helix_;

};

math::Point
get_random_point(
        util::RandomNumberGenerator &,
        int);

#endif //TEST_DOCK_MOTIF_H