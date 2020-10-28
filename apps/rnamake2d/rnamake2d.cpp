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
    //TODO add the score function cutoffs
}

void
RNAMake2D::run()  {

    static plog::ConsoleAppender< plog::TxtFormatter > console;
    plog::init(plog::info, &console);

    // TODO sometimes when NemoSampler::intialize_design() is given a difficult target it will take a long time to
    // get through it. Figure out what is going on here
    // initialize the designs
    designs_ = rnamake2d::make_designs(parameters_.dot_bracket, "",parameters_.num_designs) ;
    // find the solution for each design
    for(auto&& design : designs_) {
        // initialize the design and get initial score
        sampler_.initialize_design(design);
        design.score = vienna_fxn_.score(design);
        for(auto ii = 0; ii<parameters_.steps; ++ii) {

            if(design.score >= parameters_.score_cutoff) {
                // case to leave the loop early
                break;
            }

//            if(design.score < parameters_.sfxn_cutoff) {
//                design.score = vienna_fxn_.score(design);
//            } else {
//                design.score = sfxn_.score(design.mutant);
//            }

            sampler_.mutate(design);
            auto mutation_score(0.);

            if(design.score < parameters_.sfxn_cutoff) {
               mutation_score = vienna_fxn_.score_mutation(design);
            } else {
               mutation_score = sfxn_.score(design.mutant);
            }

            if(mc_.accept(design.score, mutation_score)) {
                design.accept(mutation_score);
            } else {
                design.reject();
            }

        }

        // probably should check if the design is above the score cutoff
        if( design.score >= parameters_.sfxn_cutoff && !parameters_.outfile.empty()) {
            // save the design information to f file
            auto outfile = std::fstream(parameters_.outfile, std::ios::app);
            outfile<<design.to_str();
            outfile.close();
        } else {
            // something close to an error message
        }
    }


}


int main( int argc, char** argv) {

    auto better_than_rnainverse = RNAMake2D();
    better_than_rnainverse.setup_options();
    CLI11_PARSE(better_than_rnainverse.app_, argc, argv);
    better_than_rnainverse.run();
}
