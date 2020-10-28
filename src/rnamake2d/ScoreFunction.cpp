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

    void
    ScoreFunction::load_weights( std::filesystem::path const& weights_file)   {
        LOGI<<"Updating score function weights... Clearing out existing weights...";
        if(!strategies_.empty())  {
            strategies_.clear();
        }

        if(!weights_.empty())  {
            weights_.clear();
        }

        if(!std::filesystem::exists(weights_file)) {
            LOGE<<"ERROR: The score function parameters file "<<weights_file.string()<<" does not exist or is inaccesible. Exiting...";
            exit(1);
        }

        for( const auto& line : base::get_lines_from_file(weights_file)) {
            const auto tokens = base::split_str_by_delimiter(line, ",");
            if(tokens.size() != 2) {
                LOGE<<"ERROR: Invalid weights file format. Exiting...";
                exit(1);
            }
            weights_.push_back(std::stof(tokens[1]));
            strategies_.emplace_back(get_strategy(tokens[0]));
        }

        LOGI<<"Loaded new parameters. Found "<<weights_.size()<<" rules.";

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

