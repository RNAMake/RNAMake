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
    //TODO add the score function cutoffs
}

void
RNAMake2D::run()  {

    static plog::ConsoleAppender< plog::TxtFormatter > console;
    plog::init(plog::info, &console);

    if(!parameters_.start_sequence.empty()) {
        if( parameters_.start_sequence.size() != parameters_.dot_bracket.size()) {
            LOGE << "ERROR: The supplied start sequence MUST be the same length as the supplied dot bracket structure. Exiting...";
            exit(1);
        } else if (parameters_.start_sequence.find_first_not_of("ACUGN") != String::npos) {
            LOGE<<"ERROR: The supplied start seqeuence can only contain the characters A, C, G, U or N. Exting...";
            exit(1);
        }
    }

    sfxn_.load_weights(parameters_.sfxn_weights);
    sfxn_.load_params(parameters_.params_dir);
    // TODO sometimes when NemoSampler::intialize_design() is given a difficult target it will take a long time to
    // get through it. Figure out what is going on here
    // initialize the designs
    // find the solution for each design
    const auto time_start = Clock::now();
    while( designs_.size() < parameters_.num_designs
            && std::chrono::duration_cast<Hours>(Clock::now() - time_start) < parameters_.max_time
    ) {
        auto design = rnamake2d::Design(parameters_.dot_bracket, parameters_.start_sequence);
        // initialize the design and get initial score
        // TODO add flag for going to score fxn
        sampler_.initialize_design(design);
        design.bp_score(vienna_fxn_.score(design));
        for(auto ii = 0; ii < parameters_.steps; ++ii) {
            // to break out of the loop, the design needs to surppass the bp_cutoff and the score_fxn cutoff
            if(design.bp_score() < parameters_.bp_cutoff) {
                // bp_cutoff not met yet
                sampler_.mutate(design);
                const auto mutation_bp = vienna_fxn_.score_mutation(design);
                if(mc_.accept(design.bp_score(), mutation_bp)) {
                    design.accept_bp(mutation_bp);
                } else {
                    design.reject();
                }
            } else {
                // the case when the bp_cutoff is already high enough
                sampler_.mutate(design);
                const auto mutation_score = sfxn_.score(design.mutant_);
                if(mc_.accept(design.score(), mutation_score)) {
                    design.accept(mutation_score);
                } else {
                    design.reject();
                }
            }
            if( design.score() >= parameters_.sfxn_cutoff && design.bp_score() >= parameters_.bp_cutoff) {
                break;
            }
        }
        // probably should check if the design is above the score cutoff
        if( design.score() >= parameters_.sfxn_cutoff
            && design.bp_score() >= parameters_.bp_cutoff) {
            // save the design information to f file
            auto outfile = std::fstream(parameters_.outfile, std::ios::app);
            outfile<<design.to_str();
            outfile.close();
            designs_.push_back(design);
        } else {
            LOGE<<"UNABLE TO GENERATE SUITABLE DESIGN. Continuing...";
            // something close to an error message
        }
    }

    if(std::chrono::duration_cast<Hours>(Clock::now() - time_start) > parameters_.max_time &&
            designs_.size() < parameters_.num_designs) {
        LOGE<<"ERROR: Had to terminate due to time constraints. Only "<<designs_.size()
                <<" designs created, not "<<parameters_.num_designs;
    }
    LOGI<<"RESULT: "<<designs_.rbegin()->to_str();
}


int main( int argc, char** argv) {

    auto better_than_rnainverse = RNAMake2D();
    better_than_rnainverse.setup_options();
    CLI11_PARSE(better_than_rnainverse.app_, argc, argv);
    better_than_rnainverse.run();
}
