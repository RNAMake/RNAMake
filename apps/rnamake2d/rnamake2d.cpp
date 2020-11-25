#include <rnamake2d/rnamake2d.h>

void
RNAMake2D::setup_options() {

    app_.add_option("target",parameters_.dot_bracket)
        ->required();

    app_.add_option("--sequence",parameters_.start_sequence);

    app_.add_option("--num_designs",parameters_.num_designs)
        ->check(CLI::PositiveNumber);

    app_.add_option("--steps",parameters_.steps)
        ->check(CLI::PositiveNumber);

    app_.add_option("--outfile", parameters_.outfile)
        ;

    app_.add_option("--max_time", parameters_.max_time)
        ->default_val(1)
        ->check(CLI::PositiveNumber);

    app_.add_option("--weights", parameters_.sfxn_weights)
        ->check(CLI::ExistingFile)
        ->default_val(base::base_dir() + "/apps/rnamake2d/settings/orig_rnamake.weights")
        ;
    app_.add_option("--params", parameters_.params_dir)
            ->check(CLI::ExistingDirectory)
            ->default_val(base::base_dir() + "/apps/rnamake2d/settings/params/")
            ;
    app_.add_flag("--cluster_output", parameters_.cluster_mode);

    //TODO add the score function cutoffs
}

void
RNAMake2D::bp_iteration_(rnamake2d::Design& design) {
    sampler_.mutate(design);
    // before doing calcs, check if this solution has been calculated already
    const auto seq = design.sequence();
    const auto is_calculated = bp_memo_.find(seq) != bp_memo_.end();
    const auto mutation_bp = is_calculated ? bp_memo_[seq] : vienna_fxn_.score_mutation(design);

    if(!is_calculated) {
        bp_memo_[seq] = mutation_bp;
    }

    if(mc_.accept(design.bp_score(), mutation_bp)) {
        design.accept_bp(mutation_bp);
    } else {
        design.reject();
    }
}

void
RNAMake2D::sfxn_iteration_(rnamake2d::Design & design) {
    sampler_.mutate(design);
    const auto seq = design.sequence();
    const auto is_calculated = memo_.find(seq) != memo_.end();
    const auto mutation_score = is_calculated ? memo_[seq] : sfxn_.score(design.mutant_);

    if(!is_calculated) {
        memo_[seq] = mutation_score;
    }

    if(mc_.accept(design.score(), mutation_score)) {
        design.accept(mutation_score);
    } else {
        design.reject();
    }
}

void
RNAMake2D::check_start_params_() const {
    if(!parameters_.start_sequence.empty()) {
        if( parameters_.start_sequence.size() != parameters_.dot_bracket.size()) {
            LOGE << "ERROR: The supplied start sequence MUST be the same length as the supplied dot bracket structure. Exiting...";
            exit(1);
        } else if (parameters_.start_sequence.find_first_not_of("ACUGN") != String::npos) {
            LOGE<<"ERROR: The supplied start seqeuence can only contain the characters A, C, G, U or N. Exting...";
            exit(1);
        }
    }
}

void
RNAMake2D::log_results_(const rnamake2d::Design & design) {
    if( design_above_threshold_(design) ) {
        // save the design information to f file
        if( parameters_.cluster_mode)  {
            std::cout<<"Target: "<<design.target()<<"\n";
            std::cout<<"Design: "<<design.sequence()<<"\n";
        } else {
            auto outfile = std::fstream(parameters_.outfile, std::ios::app);
            outfile<<design.to_str();
            outfile.close();
        }
        designs_.push_back(design);
    } else {
        LOGE<<"UNABLE TO GENERATE SUITABLE DESIGN. Continuing...";
        // something close to an error message
    }
}

void
RNAMake2D::run()  {

    static plog::ConsoleAppender< plog::TxtFormatter > console;
    if(!parameters_.cluster_mode)   {
        plog::init(plog::info, &console);
    } else {
        plog::init(plog::warning, &console);
    }

    // checking that the start params are ok
    check_start_params_();
    // loading in the parameters
    sfxn_.load_weights(parameters_.sfxn_weights);
    sfxn_.load_params(parameters_.params_dir);

    const auto time_start = timer();

    while( !exit_conditions_(time_start) ) {
        // initialize the design and get initial score
        auto design = rnamake2d::Design(parameters_.dot_bracket, parameters_.start_sequence);
        // TODO add flag for going to score fxn
        sampler_.initialize_design(design);
        design.bp_score(vienna_fxn_.score(design));
        for(auto ii = 0; ii < parameters_.steps; ++ii) {
            // to break out of the loop, the design needs to surppass the bp_cutoff and the score_fxn cutoff
            if(design.bp_score() < parameters_.bp_cutoff) {
                // bp_cutoff not met yet
                bp_iteration_(design);
            } else {
                // the case when the bp_cutoff is already high enough
                sfxn_iteration_(design);
            }

            if( design_above_threshold_(design) )  { break; }
        }
        // probably should check if the design is above the score cutoff
        log_results_(design);
    }

    exit_message_();
}


int main( int argc, char** argv) {

    auto better_than_rnainverse = RNAMake2D();
    better_than_rnainverse.setup_options();
    CLI11_PARSE(better_than_rnainverse.app_, argc, argv);
    better_than_rnainverse.run();
}
