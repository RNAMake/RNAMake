//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_SEARCH_H
#define RNAMAKE_NEW_SEARCH_H

#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_search/problem.h>

namespace motif_search {

struct Solution {
    inline
    Solution(
            motif_data_structure::MotifStateGraphOP n_graph,
            float n_score):
            graph(n_graph),
            score(n_score) {}

    motif_data_structure::MotifStateGraphOP graph;
    float score;
};

typedef std::shared_ptr<Solution> SolutionOP;


class Search {
public:

public:
    Search() {}

    ~Search() {}

    virtual
    Search *
    clone() const = 0;

public:
    virtual
    void
    setup(
            ProblemOP) = 0;


    virtual
    void
    start() = 0;

    virtual
    bool
    finished() = 0;

    virtual
    SolutionOP
    next() = 0;

};

typedef std::shared_ptr<Search> SearchOP;

}


#endif //RNAMAKE_NEW_SEARCH_H
