//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_SEARCH_H
#define RNAMAKE_NEW_SEARCH_H

#include <motif_search/problem.h>

namespace motif_search {

struct Solution {

};

template<typename SolutionType>
class Search {
public:
    typedef std::shared_ptr<SolutionType> SolutionTypeOP;

public:
    Search() {}

    ~Search() {}

public:
    virtual
    void
    setup(
            ProblemOP) = 0;


    virtual
    void
    start() = 0;

    virtual
    void
    finished() = 0;

    virtual
    SolutionTypeOP
    next() = 0;





};

}


#endif //RNAMAKE_NEW_SEARCH_H
