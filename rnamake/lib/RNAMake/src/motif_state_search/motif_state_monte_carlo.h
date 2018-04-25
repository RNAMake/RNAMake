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
        using_lookup_ = 0;
        setup_options();
    }

public:
    void
    setup(
            MotifStateGraphOP msg,
            int,
            int,
            int,
            int,
            bool);

    void
    run();

    void
    start();

    MotifGraphOP
    next();

    bool
    finished();

public:

    inline
    void
    lookup(
            StericLookup const & sl) {
        using_lookup_ = 1;
        lookup_ = sl;
    }

private:

    double
    get_score(
            BasepairStateOP);

    float
    perform_motif_swap(
            float);

    bool
    _steric_clash(
            MotifStateGraphOP);

    bool
    _seen_solution(
            MotifStateGraphOP);

protected:

    void
    setup_options();

    void
    update_var_options();

public: //option wrappers

    inline
    Options &
    options() { return options_; }

    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }

    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }

    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }

    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }

    inline
    bool
    has_option(String const & name) { return options_.has_option(name); }

    template<typename T>
    void
    set_option_value(
            String const & name,
            T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }

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
    int stage_;
    int step_;
    bool target_an_aligned_end_;
    StringIntMap seen_;
    bool enumerating_;
    int using_lookup_;
    // options
    bool sterics_;
    int max_solutions_;
    int stages_;
    int steps_;
    float accept_score_;
};


#endif //TEST_MOTIF_STATE_MONTE_CARLO_H
