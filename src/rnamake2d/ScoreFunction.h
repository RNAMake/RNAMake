#ifndef __RNAMAKE_SCORE_FUNCTION_H__
#define __RNAMAKE_SCORE_FUNCTION_H__

#include <map>
#include <filesystem>

#include <base/types.h>
#include <rnamake2d/design.h>
#include <rnamake2d/Rule.h>
#include <eternabot/strategy.h>
#include <vienna/vienna.h>

#include <base/file_io.h>
#include <base/string.h>
#include <plog/Log.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class ScoreFunction {
        //eternabot::StrategyOPs strategies_;
        std::vector<Strategy2DOP> strategies_;
        Reals weights_;
        int num_strats;

    public:
        ScoreFunction() {
            strategies_ = { } ;
            weights_ = { };
            num_strats = strategies_.size();
        }

        public:
            double
            score( Feature2DOP const & ) const ;

        public:
            void
            load_params( std::filesystem::path const& ) ;

        public:
            void
            load_weights( std::filesystem::path const& ) ;

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



#endif // __RNAMAKE_SCORE_FUNCTION_H__