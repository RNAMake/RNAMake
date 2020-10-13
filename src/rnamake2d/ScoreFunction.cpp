#include <rnamake2d/ScoreFunction.h>

namespace rnamake2d {

    double
    ScoreFunction::score(const Design & design) const {
        auto ans(0.);

        for(const auto& weight_rule : scoring_rules_) {
            ans += weight_rule.first*weight_rule.second->score(design);
        }

        return ans;

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
} // namespace rnamake2d

