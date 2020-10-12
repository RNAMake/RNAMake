#ifndef __RNAMAKE_2D_H__
#define __RNAMAKE_2D_H__

#include <base/types.h>
#include <base/application.hpp>
#include <rnamake2d/Design.h>
#include <rnamake2d/ScoreFunction.h>
#include <rnamake2d/NemoSampler.h>
#include <vienna/vienna.h>
#include <util/monte_carlo.h>


#include <CLI/CLI.hpp>
#include <nemo/nemo.cpp>

class RNAMake2D : public base::Application {
public:
    RNAMake2D() :
        app_("RNAMake2D"),
        sampler_()
    {

    }

public:
    void
    run() override;

public:
    void
    setup_options() override;

private:
    struct Parameters {
        String dot_bracket;
        int seed{1996};
        int steps{200};
        double sfxn_cutoff{80.f};
        double score_cutoff{95.f};
    };

private:
    Parameters parameters_;
    Designs designs_;
    rnamake2d::ScoreFunction sfxn_;
    rnamake2d::ViennaScore vienna_fxn_;
    rnamake2d::NemoSampler sampler_;
    util::MonteCarlo mc_;  // don't forget to seed this

public:
    CLI::App app_;
};

#endif // defined(__RNAMAKE_2D_H__)
