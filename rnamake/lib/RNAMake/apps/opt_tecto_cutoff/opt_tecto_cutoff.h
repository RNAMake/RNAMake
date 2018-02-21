//
// Created by Joseph Yesselman on 2/19/18.
//

#ifndef TEST_OPT_TECTO_CUTOFF_H
#define TEST_OPT_TECTO_CUTOFF_H

#include <stdio.h>

//RNAMake Headers
#include "base/application.hpp"


struct Construct {
    String fseq;
    String fss;
    String cseq;
    String css;
    double dg;
    int pos;
};

struct CompareConstructDG {
    inline
    bool
    operator() (
            Construct const & c1,
            Construct const & c2) const {
        return c1.dg < c2.dg;
    }

};

typedef std::vector<Construct> Constructs;

class OptTectoCutoff : public Application {
public:
    OptTectoCutoff();

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

private:
    void
    _setup();

    double
    _score(
            std::array<Real2, 6> const &);

    void
    _vary_constraints(
            std::array<Real2, 6> const &,
            std::array<Real2, 6> &,
            RandomNumberGenerator &);

    void
    _divide_dataset();

    void
    _get_scored_dataset();

    int
    _parse_constraint_position(
            String const &);

    void
    _get_initial_constraints(
            RandomNumberGenerator &);

private:
    std::vector<SixDHistogram> histos_;
    Constructs constructs_;
    std::array<Real2, 6> constraints_, new_constraints_;
    std::vector<double> exp_dgs_, norm_exp_dgs_, pred_dgs_;
    std::vector<double> avg_hit_counts_;
    double r_, avg_diff_;



};

#endif //TEST_OPT_TECTO_CUTOFF_H
