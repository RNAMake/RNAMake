//
//  motif_tree_state_path_refiner.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 3/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <sstream>
#include "motif_tree_state_path_refiner.h"

MotifTreeStateSearchSolutionOPs
MotifTreeStatePathRefiner::find_path(
    BasepairStateOP const & start,
    BasepairStateOP const & end) {
    
    // re init search
    search_ = MotifTreeStateSearch();
    set_search_options();
    MotifTreeStateSearchScorerOP scorer ( new MTSS_Astar() );
    search_.set_numeric_option("max_n_solutions", 10);
    MotifTreeStateSearchSolutionOPs solutions = search_.search(start, end, NULL);

    if(solutions.size() == 0) { return solutions; }
    std::stringstream ss;

    int i = 0;
    for (auto const & s : solutions) {
        ss << "solutions." << i << ".pdb";
        s->to_mtst().to_pdb(ss.str());
        ss.str("");
        i++;
    }
    return solutions;
    
    scorer->ss_score_weight(0);
    
    float level_weight = 1.0;
    float ss_score = 1.0;
    int count = 0;
    
    solutions.push_back(NULL);
    scorer->level_weight(2);
    while(solutions.size() > 0) {
        scorer->ss_score_weight(ss_score);
        search_.reset();
        solutions = search_.search(start, end, NULL, scorer);
        
        if(solutions.size() == 0) { break; }
        solutions[0]->to_mtst().to_pdb("solution.min.score.pdb");
        ss_score += 1.0;
        std::cout << ss_score << std::endl;
        if(ss_score > 4) {
            break;
        }
    }
    
    std::cout << level_weight << " " << ss_score << std::endl;
    
    
    
    return solutions;
}

void
MotifTreeStatePathRefiner::set_search_options() {
    for (auto const & kv : options_.numeric_options()) {
        search_.set_numeric_option(kv.first, kv.second);
    }
}