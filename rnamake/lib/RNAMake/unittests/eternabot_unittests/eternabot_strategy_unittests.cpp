//
//  eternabot_strategy_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//



#include "secondary_structure/secondary_structure_parser.h"
#include "eternabot/feature_generator.h"
#include "eternabot/strategy.h"
#include "eternabot/strategy/a_basic_test.h"
#include "eternabot/strategy/clear_plot.h"
#include "eternabot/strategy/berex_test.h"
#include "eternabot/strategy/num_of_yellow.h"

#include "eternabot_strategy_unittests.h"
#include "build/build_secondary_structure.h"
#include "build/build_motif_graph.h"
#include "secondary_structure_unittests/util.h"


namespace unittests {
namespace eternabot_unittests {

    
int
EternabotStrategyUnittest::test_creation() {
    auto strategy = eternabot::ABasicTest();
    
    return 0;
}
    
int
EternabotStrategyUnittest::test_pose_helix() {
    
    auto builder2 = BuildMotifGraph();
    auto mg = builder2.build(3);
    mg->replace_ideal_helices();
    auto dss = mg->designable_secondary_structure();
    unittests::sstruct_unittests::fill_basepairs_in_ss(dss);
    
    auto helices = dss->helices();
    if(helices.size() != 2) {
        throw UnittestException("did not get the correct number of helices");
    }
    
    return 0;
}
    
int
EternabotStrategyUnittest::test_score() {
    auto builder  = BuildSecondaryStructure();
    auto p = builder.build_helix(10);
    auto generator = eternabot::FeatureGenerator();
    
    auto features = generator.get_features(p);
    
    auto strategy1 = eternabot::ABasicTest();
    auto strategy2 = eternabot::CleanPlotStackCapsandSafeGC();
    auto strategy3 = eternabot::BerexTest();
    auto strategy4 = eternabot::NumofYellowNucleotidesperLengthofString();
    
    auto score1 = strategy1.score(features);
    auto score2 = strategy2.score(features);
    auto score3 = strategy3.score(features);
    auto score4 = strategy4.score(features);
    
    //std::cout << score1 << " " << score2 << " " << score3 << " " << score4 << std::endl;
    
    return 0;
}
    
int
EternabotStrategyUnittest::test_score_compare() {
    auto seq = "CUCGAUAAACUAAGCUGUCCAAAAGGACAGCUUAGUUUAUCGAG";
    auto ss  = "((((((((((((((((((((....))))))))))))))))))))";
    
    auto parser = sstruct::SecondaryStructureParser();
    auto p = parser.parse_to_pose(seq, ss);
    
    auto generator = eternabot::FeatureGenerator();
    auto features = generator.get_features(p);
    
    auto strategy1 = eternabot::ABasicTest();
    auto strategy2 = eternabot::CleanPlotStackCapsandSafeGC();
    auto strategy3 = eternabot::BerexTest();
    auto strategy4 = eternabot::NumofYellowNucleotidesperLengthofString();
    
    auto score1 = strategy1.score(features);
    auto score2 = strategy2.score(features);
    auto score3 = strategy3.score(features);
    auto score4 = strategy4.score(features);
    
    std::cout << score1 << " " << score2 << " " << score3 << " " << score4 << std::endl;
    
    
    return 0;
}

    
int
EternabotStrategyUnittest::run() {
    test_creation();
    test_pose_helix();
    //test_score();
    test_score_compare();
    return 0;
}
    
int
EternabotStrategyUnittest::run_all() {
    return 0;
}

    
    
}
}