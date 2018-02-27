//
// Created by Joseph Yesselman on 11/4/17.
//

#ifndef TEST_MOTIF_STATE_MONTE_CARLO_H
#define TEST_MOTIF_STATE_MONTE_CARLO_H

#include <base/option.h>
#include <util/steric_lookup.hpp>
#include <util/monte_carlo.h>
#include <motif_data_structures/motif_state_tree.h>
#include <motif_data_structures/motif_state_graph.hpp>

class MotifStateMonteCarlo {
public:
    MotifStateMonteCarlo(
            std::vector<MotifStateOPs> const & mses):
            mses_(mses),
            mc_(MonteCarlo(0.5f)),
            rng_(RandomNumberGenerator()),
            options_(Options()){
    }

public:
    void
    setup(
            MotifStateGraphOP msg,
            int,
            int,
            int,
            int);

    void
    run();

private:

    double
    get_score(
            BasepairStateOP);

    void
    perform_motif_swap();


private:
    MonteCarlo mc_;
    RandomNumberGenerator rng_;
    StericLookup lookup_;
    Points beads_;
    Options options_;
    std::vector<MotifStateOPs> mses_;
    BasepairStateOP end_, end_flip_;
    MotifStateOP start_m_;
    MotifStateGraphOP msg_;
    double score_, r_diff_, r_diff_flip_;
    int ni_, nj_, ei_, ej_;
    int org_num_;

};


#endif //TEST_MOTIF_STATE_MONTE_CARLO_H
