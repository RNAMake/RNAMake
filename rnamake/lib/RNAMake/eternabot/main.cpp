//
//  main.cpp
//  eternabot
//
//  Created by Joseph Yesselman on 2/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "secondary_structure_tree.h"
#include "vienna.h"
#include "sequence_designer.h"
#include "eternabot_scorer.h"
#include "cartesian_product.h"


int
test_fold() {
    
    Vienna v;
    FoldResult fr = v.vfold("GGGAAACCC");
    //std::cout << fr.structure << std::endl;
    FoldResult fr2 = v.vfold("GGGGGGGGGGGGGGGCCCCCCCCCCCCC");
    //std::cout << fr2.free_energy << std::endl;
    return 1;
}

int
test_cofold() {
    Vienna v;
    FoldResult fr = v.vcofold("GGGAAACCCCCC&AAAAAGGGAAA");
    //std::cout << fr.structure << std::endl;
    return 1;
}

int
test_design_sequence() {
    SequenceDesigner designer;
    designer.design("NNNNNUUCGNNNNN", "(((((....)))))");
    return 1;
}

int
test_design_sequence_2() {
    SequenceDesigner designer;
    String seq = designer.design("NNNNNNNNNNNNNNNNNNNNNNNNN&NNNNNNNNNNNNNNNNNNNNNNNNN",
                                 "(((((((((((((((((((((((((&)))))))))))))))))))))))))");
    std::cout << seq << std::endl;
    return 1;
}

int
test_eternabot_scorer() {
    SecondaryStructureTree sstree("(((....)))", "GUCUUCGGAC");
    EternabotScorer scorer;
    scorer.setup("(((....)))");
    for(int i = 0; i < 100000; i++) {
        if (i % 1000 == 0 ) { std::cout << i << std::endl; }
        std::cout << scorer.score_sstree(sstree) << std::endl;
        exit(0);
    }
    return 1;
}

int
enumerate_bps_in_design(String ss, String seq) {
    SecondaryStructureTree sstree(ss, seq);
    EternabotScorer scorer;
    scorer.setup(ss);
    
    Strings spl = split_str_by_delimiter("GC CG AU UA", " ");
    SecondaryStructureNodeOPs designable_bps = sstree.get_designable_bps();
    std::vector<Strings> bps;
    int i = 0;
    for(i = 0; i < designable_bps.size(); i++) { bps.push_back(spl); }
    CartesianProduct<String> combos(bps);
    Strings current;
    float best = -1000;
    String best_sequence;
    float score = 0;
    int j = 0;
    while(!combos.end()) {
        j++;
        current = combos.next();
        for(i = 0; i < designable_bps.size(); i++) {
            designable_bps[i]->bp_type(current[i]);
        }
        score = scorer.score_sstree(sstree);
        if(score > best) {
            best = score;
            SSandSeqOP seq_and_seq = sstree.get_ss_and_seq();
            best_sequence = seq_and_seq->seq;
        }
    }
    std::cout << best_sequence << " " << best << std::endl;

    return 0;
}

int
test_eternabot_enumerate_bps() {
    enumerate_bps_in_design("((((((((....))))))))", "NNNNNNNNUUCGNNNNNNNN");
    return 1;
}

int
test_sequencer_designer_2() {
    SequenceDesigner designer;
    String seq = designer.design("NNNGAUAUGGNNNNNNNNNNNNCGCAAAUNNNNNNNNNNNCGACGAAANNNNNNNNGGAAACNNNNNNNNUGGAGNNNNNNNNNNNACUCGUACGNNNNNNNNNNNNCCUAAGUCNNN",
                        "((((((..(((((((((((((((.....(((((((((((((......((((((((((....))))))))))...))))))))))))).......)))))))))))))))...))))))");
    std::cout << seq << std::endl;

    return 1;
}




int main(int argc, const char * argv[]) {
    
    if (test_fold() == 0)               { std::cout << "test_fold failed" << std::endl;  }
    if (test_cofold() == 0)             { std::cout << "test_cofold failed" << std::endl;  }
    //if (test_design_sequence() == 0)  { std::cout << "test_design_sequence failed" << std::endl;  }
    //if (test_design_sequence_2() == 0)  { std::cout << "test_design_sequence failed" << std::endl;  }
    //if (test_eternabot_scorer() == 0)   { std::cout << "test_eternabot_scorer failed" << std::endl;  }
    //test_eternabot_scorer();
    //test_eternabot_enumerate_bps();
    test_sequencer_designer_2();
    
    return 0;
}
