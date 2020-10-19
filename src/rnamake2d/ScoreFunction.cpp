#include <rnamake2d/ScoreFunction.h>

namespace rnamake2d {

    double
    ScoreFunction::score( Feature2DOP const &  feature ) const {
        auto score(0.);
        for(auto ii = 0; ii < num_strats; ++ii) {
            score += strategies_[ii]->score(feature)*weights_[ii];
        }
        return score;
    }



    double
    ViennaScore::score(Design const& Design) {

        vienna_.fold(Design.sequence);
        const auto prediction = vienna_.get_structure();

        auto ans(0);
        for( auto pred = prediction.begin(), targ = Design.target.begin();
             pred != prediction.end();
             ++pred, ++targ)  {

            if(*pred =='.' && *pred == *targ) {
                ++ans;
            } else if (*pred != '.' && *targ != '.') {
                ++ans;
            }
        }

        return 100.0*double(ans)/double(prediction.size());
    }

    double
    ViennaScore::score_mutation(const Design & Design) {

        vienna_.fold(Design.candiate);
        const auto prediction = vienna_.get_structure();

        auto ans(0);
        for( auto pred = prediction.begin(), targ = Design.target.begin();
             pred != prediction.end();
             ++pred, ++targ)  {

            if(*pred =='.' && *pred == *targ) {
                ++ans;
            } else if (*pred != '.' && *targ != '.') {
                ++ans;
            }
        }

        return 100.0*double(ans)/double(prediction.size());
    }
} // namespace rnamake2d

