//
//  scorer_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "scorer_unittest.h"

#include "eternabot/scorer.h"
#include "build/build_secondary_structure.h"
#include "secondary_structure/secondary_structure_parser.h"

namespace unittests {
namespace eternabot_unittests {
    
int
ScorerUnittest::test_creation() {
    auto scorer = eternabot::Scorer();
    auto builder = BuildSecondaryStructure();
    auto p = builder.build_hairpin(20);
    
    scorer.setup(p);
    auto score = scorer.score_secondary_structure(p);
    std::cout << p->sequence() << std::endl;
    std::cout << p->dot_bracket() << std::endl;
    
    std::cout << score << std::endl;
    
    return 0;
}
    
int
ScorerUnittest::test_score_compare() {
    auto seq = "CGCUACUUUCGACGCUGCGCAAAAGCGCAGCGUCGAAAGUAGCG";
    auto ss  = "((((((((((((((((((((....))))))))))))))))))))";
        
    auto parser = sstruct::SecondaryStructureParser();
    auto p = parser.parse_to_pose(seq, ss);
    auto scorer = eternabot::Scorer();
    
    scorer.setup(p);
    auto score = scorer.score_secondary_structure(p);
    std::cout << score << std::endl;

    
    return 0;
}
    
int
ScorerUnittest::run() {
    //test_creation();
    test_score_compare();
    return 0;
}

int
ScorerUnittest::run_all() {
    return 0;
}
    

}
}