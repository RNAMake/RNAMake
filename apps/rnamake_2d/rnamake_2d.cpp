#include <rnamake_2d/rnamake_2d.h>

void
RNAMake2D::setup_options() {
#if 0
//    app_.add_option("--dot_bracket",parameters_.dot_bracket)
//        ->required();
#endif  // 0
}

void
RNAMake2D::run()  {
    designs_.emplace_back(String{"(((...)))"});
    for(auto&& design : designs_) {
        // initialize the design
        sampler_.initialize_design(design);
        std::cout<<designs_[0].sequence<<std::endl;
        for(auto ii = 0; ii<parameters_.steps; ++ii) {

            if(design.score >= parameters_.score_cutoff) {
                // case to leave the loop early
                break;
            }

            auto candidate_score(0.);
            if(design.score < parameters_.sfxn_cutoff) {
               candidate_score = vienna_fxn_.score(design);
            } else {
               candidate_score = sfxn_.score(design);
            }

            if(mc_.accept(design.score, candidate_score)) {
                design.accept();
            } else {
                design.reject();
            }

        }

        // probably should check if the design is above the score cutoff
        if( design.score >= parameters_.sfxn_cutoff) {
            // save the design information to f file
        } else {
            // something close to an error message
        }
    }


}

int rnamake2d::Design::design_ct = 0; 

int main( int argc, char** argv) {

    auto better_than_rnainverse = RNAMake2D();
    better_than_rnainverse.setup_options();
    CLI11_PARSE(better_than_rnainverse.app_,argc,argv);
    better_than_rnainverse.run();
}
