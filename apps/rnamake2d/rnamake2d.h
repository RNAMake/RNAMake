#ifndef __RNAMAKE_2D_H__
#define __RNAMAKE_2D_H__

#include <chrono>

using Hours = std::chrono::hours;
using Clock = std::chrono::system_clock;

#include <base/types.h>
#include <base/application.hpp>
#include <base/settings.h>
#include <rnamake2d/design.h>
#include <rnamake2d/ScoreFunction.h>
#include <rnamake2d/NemoSampler.h>
#include <rnamake2d/ss_motif.h>

#include <vienna/vienna.h>
#include <util/monte_carlo.h>

#include <CLI/CLI.hpp>
#include <plog/Log.h>
#include <plog/Appenders/ConsoleAppender.h>

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
        String outfile;
        String dot_bracket;
        String start_sequence;
        int num_designs = 5;
        int seed{1996};
        int steps{200};
        double sfxn_cutoff{80.f};
        double bp_cutoff{90.f};
        Hours max_time;
        std::filesystem::path sfxn_weights;
        std::filesystem::path params_dir;
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
