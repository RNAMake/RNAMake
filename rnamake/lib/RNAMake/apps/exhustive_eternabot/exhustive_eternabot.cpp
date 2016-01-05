//
//  exhustive_eternabot.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/5/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include <fstream>
#include "exhustive_eternabot.h"

#include "util/cartesian_product.h"
#include "eternabot/scorer.h"
#include "secondary_structure/secondary_structure_parser.h"



int main(int argc, const char * argv[]) {
    auto seq = "GGGGAAACCCCGAGC";
    auto ss  = "((((...))))....";
    
    auto parser = sstruct::SecondaryStructureParser();
    auto p = parser.parse_to_pose(seq, ss);
    
    auto pairs = std::vector<Strings>();
    pairs.push_back(Strings{"A", "U"});
    pairs.push_back(Strings{"U", "A"});
    pairs.push_back(Strings{"C", "G"});
    pairs.push_back(Strings{"G", "C"});
    //pairs.push_back(Strings{"G", "U"});
    //pairs.push_back(Strings{"U", "G"});

    
    auto all_pairs = std::vector<std::vector<Strings>>(p->basepairs().size());
    for(int i = 0; i < all_pairs.size(); i++) {
        all_pairs[i] = pairs;
    }

    auto pair_iterator = CartesianProduct<Strings>(all_pairs);
    auto current = std::vector<Strings>();
    auto bps = p->basepairs();

    auto scorer = eternabot::Scorer();
    scorer.setup(p);
    float score = 0;
    float best = 0;
    auto best_seq = p->sequence();
    
    auto out = std::ofstream("4bp_cccc.out");
    
    int i = 0;
    while (!pair_iterator.end()) {
        current = pair_iterator.next();
        i = 0;
        for(auto const & p : current) {
            bps[i]->res1()->name(p[0]);
            bps[i]->res2()->name(p[1]);
            i++;
        }
        score =  scorer.score_secondary_structure(p);
        if(score > best) {
            best = score;
            best_seq = p->sequence();
        }
        //std::cout << p->sequence() << " " << scorer.score_secondary_structure(p) << std::endl;
        out << p->sequence() << " " << score << std::endl;
        
    }
    
    std::cout << best << " " << best_seq << std::endl;
    out.close();
    
    return 0;
}


