//
//  scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__scorer__
#define __RNAMake__scorer__

#include <stdio.h>

#include "eternabot/feature_generator.h"
#include "eternabot/strategy.h"
#include "eternabot/strategy/a_basic_test.h"
#include "eternabot/strategy/clear_plot.h"
#include "eternabot/strategy/berex_test.h"
#include "eternabot/strategy/num_of_yellow.h"
#include "eternabot/strategy/direction_of_gc.h"
#include "eternabot/strategy/modified_berex_test.h"

namespace eternabot {

class StrategyFactory {
public:
    StrategyOP
    get_strategy(
            String const & name) {
        if(name == "ABasicTest") {
            return std::make_shared<ABasicTest>();
        }
        else if(name == "CleanPlotStackCapsandSafeGC") {
            return std::make_shared<CleanPlotStackCapsandSafeGC>();
        }
        else if(name == "DirectionofGCPairsinMultiLoops") {
            return std::make_shared<DirectionofGCPairsinMultiLoops>();
        }
        else if(name == "BerexTest") {
            return std::make_shared<BerexTest>();
        }
        else if(name == "NumofYellowNucleotidesperLengthofString") {
            return std::make_shared<NumofYellowNucleotidesperLengthofString>();
        }
        else if(name == "ModifiedBerexTest") {
            return std::make_shared<ModifiedBerexTest>();
        }
        else {
            throw secondary_structure::Exception("unknown strategy name");
        }
    }
};

class Scorer {
public:
    
    Scorer() :
        generator_( FeatureGenerator() ),
        strategies_( StrategyOPs() ),
        weights_ ( Floats() ) {

        auto strat_factory = StrategyFactory();
        auto names = Strings{"ABasicTest", "CleanPlotStackCapsandSafeGC", "DirectionofGCPairsinMultiLoops", "ModifiedBerexTest",
                             "NumofYellowNucleotidesperLengthofString"};
        weights_ = Floats{0.09281782, 0.1250677, 0.2156337, 0.3661276, 0.2230357};
        for(auto const & name : names) {
            strategies_.push_back(strat_factory.get_strategy(name));
        }

        scores_ = Floats( strategies_.size() );
        
    }

    Scorer(
            Strings const & strategy_names,
            Floats const & weights):
            weights_(weights) {
        auto strat_factory = StrategyFactory();
        for(auto const & name : strategy_names) {
            strategies_.push_back(strat_factory.get_strategy(name));
        }

        scores_ = Floats( strategies_.size() );
    }


    ~Scorer() {}

public:

    void
    setup(secondary_structure::PoseOP const &);
    
    float
    score_secondary_structure(secondary_structure::PoseOP const &);

    float
    print_scores(secondary_structure::PoseOP const &);

    Floats
    get_scores(secondary_structure::PoseOP const &);

    
public:
    
    Floats const &
    scores() { return scores_; }

    Strings
    strategy_names() {
        auto names = Strings();
        for(auto const & s : strategies_) {
            names.push_back(s->name());
        }
        return names;
    }

    FeaturesOP
    features() { return features_; }

private:
    FeatureGenerator generator_;
    FeaturesOP features_;
    StrategyOPs strategies_;
    Floats weights_, scores_;
    float mean_, stdev_, total_score_ = 0.0;

};

typedef std::shared_ptr<Scorer> ScorerOP;

//ScorerOP
//classic_scorer();


}

#endif /* defined(__RNAMake__scorer__) */
