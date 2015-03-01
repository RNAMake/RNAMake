//
//  motif_tree_state_search_solution.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_search_solution.h"
#include "motif_tree_state.h"
#include "motif.h"
#include "xyzMatrix.h"

MotifTreeStateSearchSolution::MotifTreeStateSearchSolution(
    MotifTreeStateSearchNodeOP const & node,
    float score):
    score_ (score )
{
    path_ = MotifTreeStateSearchNodeOPs();
    MotifTreeStateSearchNodeOP current = node;
    while ( current != NULL) {
        path_.push_back(current);
        current = current->parent();
    }
    
    std::reverse(path_.begin(), path_.end());
    
}


MotifTreeStateTree
MotifTreeStateSearchSolution::to_mtst() {
    MotifTreeStateOP start_mts = path_[0]->mts();
    Motif rmotif = ref_motif();
    Matrix r;
    dot(start_mts->end_states()[0]->r_T(), rmotif.ends()[0]->r(), r);
    Point trans = -rmotif.ends()[0]->d();
    Transform t (r, trans);
    rmotif.transform(t);
    Point bp_pos_diff = start_mts->end_states()[0]->d() - rmotif.ends()[0]->d();
    rmotif.move(bp_pos_diff);
    MotifTreeStateOP final_start_mts (new MotifTreeState("mtss_start", start_mts->start_index(), start_mts->size(), start_mts->score(), start_mts->beads(), start_mts->end_states(), start_mts->flip(), rmotif.to_str()) );
    MotifTreeStateTree mtst;
    mtst.sterics(0);
    int i = -1;
    int success = 0, same = 1;
    float diff = 0;
    for( auto const & n : path_) {
        i++;
        if(i == 0) { continue; }
        success = 0;
        for(int j = 0; j < mtst.last_node()->states().size(); j++) {
            MotifTreeStateNodeOP node = mtst.add_state(n->mts(), NULL, j);
            if(node == NULL) { continue; }
            same = 1;
            for (int k = 0; k < node->states().size(); k++) {
                if(n->active()[k] == 0 && node->states()[k] == NULL) { continue; }
                if(n->active()[k] == 0 && node->states()[k] != NULL) {
                    same = 0; break;
                }
                if(n->active()[k] == 1 && node->states()[k] == NULL) {
                    same = 0; break;
                }
                diff = n->states()[k]->diff(node->states()[k]);
                if(diff > 0.1) {
                    std::cout << "made it" << std::endl;
                    same = 0; break;
                }
            }
            if(same) {
                success = 1;
                break;
            }
            else {
                mtst.remove_node(node);
            }
        }
        
        if(!success) {
            std::cout << "fail" << std::endl;
            exit(0);
        }
    }
    return mtst;
}
