//
//  eternabot_strategy_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//



#include "eternabot/feature_generator.h"
#include "eternabot/strategy.h"
#include "eternabot/strategy/a_basic_test.h"

#include "eternabot_strategy_unittests.h"
#include "build/build_secondary_structure.h"


namespace unittests {
namespace eternabot {

    
int
EternabotStrategyUnittest::test_creation() {
    auto strategy = ::eternabot::ABasicTest();
    
    return 0;
}
    

int
EternabotStrategyUnittest::test_score() {
    auto builder  = BuildSecondaryStructure();
    auto p = builder.build_helix(10);
    auto generator = ::eternabot::FeatureGenerator();
    
    auto features = generator.get_features(p);
    
    auto strategy1 = ::eternabot::ABasicTest();
    
    auto score = strategy1.score(features);
    
    
    return 0;
}
    
    
int
EternabotStrategyUnittest::run() {
    test_creation();
    test_score();
    return 0;
}
    
int
EternabotStrategyUnittest::run_all() {
    return 0;
}

    
    
}
}