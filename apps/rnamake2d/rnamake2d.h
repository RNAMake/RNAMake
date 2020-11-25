#ifndef __RNAMAKE_2D_H__
#define __RNAMAKE_2D_H__

#include <chrono>

using Hours = std::chrono::hours;
using Clock = std::chrono::system_clock;
constexpr auto timer = Clock::now;

#include <base/types.h>
#include <base/application.hpp>
#include <base/settings.h>
#include <rnamake2d/design.h>
#include <rnamake2d/score_function.h>
#include <rnamake2d/nemo_sampler.h>
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
    bool
    design_above_threshold_(const rnamake2d::Design& design ) const {
        return design.score() >= parameters_.sfxn_cutoff && design.bp_score() >= parameters_.bp_cutoff;
    }
private:
    void
    log_results_(const rnamake2d::Design &);

private:
    void
    exit_message_() const  {
        if(std::chrono::duration_cast<Hours>(Clock::now() - timer()) > parameters_.max_time &&
           designs_.size() < parameters_.num_designs) {
            LOGE<<"ERROR: Had to terminate due to time constraints. Only "<<designs_.size()
                <<" designs created, not "<<parameters_.num_designs;
        }
        LOGI<<"RESULT: "<<designs_.rbegin()->to_str();
    }

private:
    bool
    exit_conditions_(const decltype(timer()) & time_start) const {
        return std::chrono::duration_cast<Hours>(timer() - time_start) < parameters_.max_time &&
            designs_.size() >= parameters_.num_designs ;
    }

private:
    void
    check_start_params_() const ;

private:
    void
    bp_iteration_(rnamake2d::Design&);

private:
    void
    sfxn_iteration_(rnamake2d::Design&);

private:
    struct Parameters {
        String outfile;
        String dot_bracket;
        String start_sequence;
        int num_designs = 1;
        int seed{1996};
        int steps{2000};
        double sfxn_cutoff{80.f};
        double bp_cutoff{95.f};
        Hours max_time;
        std::filesystem::path sfxn_weights;
        std::filesystem::path params_dir;
        bool cluster_mode = false;
    };

private:
    Parameters parameters_;
    Designs designs_;
    rnamake2d::ScoreFunction sfxn_;
    rnamake2d::ViennaScore vienna_fxn_;
    rnamake2d::NemoSampler sampler_;
    util::MonteCarlo mc_;  // don't forget to seed this
    std::unordered_map<String, float> bp_memo_;
    std::unordered_map<String, float> memo_;

public:
    CLI::App app_;
};

#endif // defined(__RNAMAKE_2D_H__)
