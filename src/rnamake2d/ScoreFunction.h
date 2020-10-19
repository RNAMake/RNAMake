#ifndef __RNAMAKE_SCORE_FUNCTION_H__
#define __RNAMAKE_SCORE_FUNCTION_H__

#include <map>

#include <base/types.h>
#include <rnamake2d/Design.h>
#include <rnamake2d/Rule.h>
#include <eternabot/strategy.h>
#include <vienna/vienna.h>

// actual strategies...

// first those from eternabot
#include <eternabot/strategy/berex_test.h>
#include <eternabot/strategy/clear_plot.h>
#include <eternabot/strategy/direction_of_gc.h>
#include <eternabot/strategy/modified_a_basic_test.h>
#include <eternabot/strategy/modified_berex_test.h>
#include <eternabot/strategy/modified_clear_plot.h>
#include <eternabot/strategy/modified_direction_of_gc.h>
#include <eternabot/strategy/modified_num_of_yellow.h>
#include <eternabot/strategy/num_of_yellow.h>

namespace rnamake2d {

    class ScoreFunction {
        eternabot::StrategyOPs strategies_;
        Reals weights_;
        int num_strats;
    public:
        ScoreFunction() {

            strategies_ = {
                //std::make_shared<eternabot::BerexTest>(),
                std::make_shared<eternabot::CleanPlotStackCapsandSafeGC>(),
                std::make_shared<eternabot::DirectionofGCPairsinMultiLoops>(),
                std::make_shared<eternabot::ModifiedABasicTest>(),
                std::make_shared<eternabot::ModifiedBerexTest>(),
                //std::make_shared<eternabot::ModifiedNumofYellowNucleotidesperLengthofString>(),
                std::make_shared<eternabot::NumofYellowNucleotidesperLengthofString>()
            } ;
            // weights
            weights_ = {
                    0.1250677,
                    0.2156337,
                    0.09281782,
                    0.3661276,
                    0.2230357
            };

            num_strats = strategies_.size();
        }

        public:
            double
            score( Feature2DOP const & ) const ;


        private:
            std::map<double,RuleOP> scoring_rules_;
    };

    class ViennaScore {
        private:
            vienna::Vienna vienna_;
        public:
            double
            score(Design const&);

            double
            score_mutation(Design const&);
    };

}



#endif // __RNAMAKE_SCORE_FUNCTION_H__src/eternabot/strategy/a_basic_test.h
