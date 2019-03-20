//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_PROBLEM_H
#define RNAMAKE_NEW_PROBLEM_H

#include <memory>

//RNAMake includes
#include <util/steric_lookup.hpp>
#include <structure/basepair_state.h>

namespace motif_search {

struct Problem {
    inline
    Problem(
            structure::BasepairStateOP n_start,
            structure::BasepairStateOP n_end,
            util::StericLookupNewOP n_lookup,
            bool n_target_an_aligned_end):
            start(n_start),
            end(n_end),
            lookup(n_lookup),
            target_an_aligned_end(n_target_an_aligned_end) {}

    structure::BasepairStateOP start, end;
    util::StericLookupNewOP lookup;
    bool target_an_aligned_end;
};

typedef std::shared_ptr<Problem> ProblemOP;

}


#endif //RNAMAKE_NEW_PROBLEM_H
