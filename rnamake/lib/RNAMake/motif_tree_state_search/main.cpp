//
//  main.cpp
//  motif_tree_state_search
//
//  Created by Joseph Yesselman on 2/21/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <random>
#include "motif_tree_state_selector.h"
#include "motif_tree_state_search_scorer.h"
#include "motif_tree_state_search_node.h"
#include "motif_tree_state_library.h"
#include "motif_tree_state_tree.h"
#include "motif_tree_state_search.h"

MotifTreeStateTree
get_two_mts_tree(int size=2) {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    MotifTreeStateTree mtst;
    
    srand(unsigned(time(NULL)));
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0,1);
    
    while ( mtst.nodes().size() < size + 1) {
        int pos = dist(mt)*mts_lib.motif_tree_states().size();
        MotifTreeStateOP mts = mts_lib.motif_tree_states()[pos];
        mtst.add_state(mts, NULL);
    }
    return mtst;
}

MotifTreeStateTree
get_two_way_helix_mts_tree(int size=10) {;
    std::vector<MotifTreeStateLibrary> mts_libs(2);
    mts_libs[0] = MotifTreeStateLibrary ( HELIX );
    mts_libs[1] = MotifTreeStateLibrary ( TWOWAY );
    MotifTreeStateTree mtst;
    
    srand(unsigned(time(NULL)));
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0,1);

    int i = 0;
    int pos = 0;
    int count = 0;
    while (i < size) {
        if( i % 2 == 0) { pos = 0; }
        else            { pos = 1; }
        
        int mts_pos = dist(mt)*(mts_libs[pos].motif_tree_states().size()-1);
        MotifTreeStateOP mts = mts_libs[pos].motif_tree_states()[mts_pos];
        
        MotifTreeStateNodeOP node = mtst.add_state(mts, NULL);
        if (node != NULL) { i++; }
        count++;
        if(count > 10000) { break; }
        
    }
     
    return mtst;
    
}

int
test_selector() {
    MotifTreeStateSelector selector = default_selector(MotifTypes());
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    MotifTreeStateSearchNodeOP node ( new MotifTreeStateSearchNode( mts_lib.motif_tree_states()[0], NULL, -1));
    MTSTypePairs children = selector.get_children_mts(node);
    return 1;
}

int
test_scorer() {
    MotifTreeStateTree mtst = get_two_mts_tree();
    BasepairStateOP target = mtst.last_node()->active_states()[0];
    MTSS_GreedyBestFirstSearch scorer;
    scorer.setup(target);
    BasepairStateOP start = mtst.nodes()[0]->states()[0];
    MotifTreeStateOP mts ( new MotifTreeState( ref_mts() ));
    MotifTreeStateSearchNodeOP node ( new MotifTreeStateSearchNode(mts, NULL, -1));
    return 1;
}

int
test_search() {
    MotifTreeStateSelectorOP selector (new MotifTreeStateSelector(default_selector(MotifTypes(), "all")));
    MotifTreeStateSearch mtss;
    mtss.set_numeric_option("accept_score", 1.0);
    mtss.set_numeric_option("max_node_level", 2);
    
    for (int i = 0; i < 1; i++) {
        MotifTreeStateTree mtst = get_two_mts_tree();
        mtss.reset();
        MotifTreeStateSearchSolutionOPs solutions = mtss.search(mtst.nodes()[0]->states()[0], mtst.nodes().back()->active_states()[0], selector);
        
        /*for(int j = 1; j < 3; j++) {
            std::cout << mtst.nodes()[j]->mts()->name() << " " << solutions[0]->path()[j]->mts()->name() << std::endl;
        }*/
        
        if(solutions.size() == 0) {
            std::cout << solutions.size() << std::endl;
         
        }
    }
    return 1;
}

int
test_search2() {
    MotifTreeStateSelectorOP selector (new MotifTreeStateSelector(default_selector(MotifTypes(), "all")));
    MotifTreeStateSearch mtss;
    mtss.set_numeric_option("accept_score", 10.0);
    mtss.set_numeric_option("max_node_level", 10);
    mtss.set_numeric_option("max_n_solutions", 1);
    
    for (int i = 0; i < 100; i++) {
        MotifTreeStateTree mtst = get_two_mts_tree(10);
        //mtst.to_pdb("start.pdb");
        mtss.reset();
        MotifTreeStateSearchSolutionOPs solutions = mtss.search(mtst.nodes()[0]->states()[0], mtst.nodes().back()->active_states()[0], selector);

        std::cout << solutions.size() << std::endl;
        //MotifTreeStateTree mtst2 = solutions[0]->to_mtst();
        //mtst2.to_pdb("solution.pdb");
    }
    return 1;
}

int
test_search3() {
    MotifTreeStateSelectorOP selector (new MotifTreeStateSelector(default_selector(MotifTypes())));
    MotifTreeStateSearchScorerOP scorer ( new MTSS_Astar() );
    MotifTreeStateSearchScorerOP scorer2 ( new MTSS_GreedyBestFirstSearch() );

    MotifTreeStateSearch mtss;
    mtss.set_numeric_option("accept_score", 10.0);
    mtss.set_numeric_option("max_node_level", 100);
    mtss.set_numeric_option("max_n_solutions", 1);
    
    for (int i = 0; i < 1; i++) {
        MotifTreeStateTree mtst = get_two_way_helix_mts_tree(50);
        mtst.to_pdb("start.pdb");
        mtss.reset();
        MotifTreeStateSearchSolutionOPs solutions = mtss.search(mtst.nodes()[0]->states()[0], mtst.nodes().back()->active_states()[0], selector, scorer);
        
        std::cout << solutions.size() << std::endl;
        if(solutions.size() > 0) {
            MotifTreeStateTree mtst2 = solutions[0]->to_mtst();
            mtst2.to_pdb("solution.pdb");
        }
        else {
            solutions = mtss.search(mtst.nodes()[0]->states()[0], mtst.nodes().back()->active_states()[0], selector, scorer2);
            std::cout << solutions.size() << std::endl;
            MotifTreeStateTree mtst2 = solutions[0]->to_mtst();
            mtst2.to_pdb("solution.pdb");
        }
    }
    return 1;
}



int main(int argc, const char * argv[]) {
    if (test_selector() == 0)     { std::cout << "test_selector failed" << std::endl;  }
    if (test_scorer() == 0)       { std::cout << "test_scorer failed" << std::endl;  }
    if (test_search3() == 0)       { std::cout << "test_search failed" << std::endl;  }

    return 0;
}











