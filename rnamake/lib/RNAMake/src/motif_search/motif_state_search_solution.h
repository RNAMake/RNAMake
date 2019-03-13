//
//  motif_state_search_solution.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_search_solution__
#define __RNAMake__motif_state_search_solution__

#include <stdio.h>

//RNAMake Headers
#include "motif_data_structure/motif_state_tree.h"
#include "motif_search/motif_state_search_node.h"

namespace motif_search {

class MotifStateSearchSolution {
public:
    MotifStateSearchSolution(
            MotifStateSearchNodeOP const & node,
            float score) :
            score_(score) {
        _get_path(node);
    }

    ~MotifStateSearchSolution() {}

public:

    motif_data_structure::MotifStateTreeOP
    to_mst();

    motif_data_structure::MotifTreeOP
    to_motif_tree() {
        return to_mst()->to_motif_tree();
    }

    void
    to_pdb(
            String const & fname,
            int renumber) {
        return to_motif_tree()->to_pdb(fname, renumber);
    }

public:
    inline
    float
    score() { return score_; }

    inline
    MotifStateSearchNodeOPs const &
    path() { return path_; }

private:

    void
    _get_path(
            MotifStateSearchNodeOP const &);

private:
    MotifStateSearchNodeOPs path_;
    float score_;

};

typedef std::shared_ptr<MotifStateSearchSolution> MotifStateSearchSolutionOP;
typedef std::vector<MotifStateSearchSolutionOP> MotifStateSearchSolutionOPs;

}

#endif /* defined(__RNAMake__motif_state_search_solution__) */
