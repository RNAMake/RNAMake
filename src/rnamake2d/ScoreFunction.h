#ifndef __RNAMAKE_SCORE_FUNCTION_H__
#define __RNAMAKE_SCORE_FUNCTION_H__

#include <map>

#include <base/types.h>
#include <rnamake2d/Design.h>
#include <rnamake2d/Rule.h>

#include <vienna/vienna.h>

namespace rnamake2d {

    class ScoreFunction {

        public:
            double
            score( Design const& ) const ;

            double
            score_mutation( Design const& ) const ;

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