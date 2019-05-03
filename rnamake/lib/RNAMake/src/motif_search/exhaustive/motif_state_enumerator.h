//
// Created by Joseph Yesselman on 2019-03-30.
//

#ifndef RNAMAKE_NEW_MOTIF_STATE_ENUMERATOR_H
#define RNAMAKE_NEW_MOTIF_STATE_ENUMERATOR_H

#include <motif_search/solution_topology.h>
#include <motif/motif_state_aligner.h>

namespace motif_search {
namespace exhaustive {

/*class NoStericAligner : motif::MotifStateAligner {

};
 */

class MotifStateEnumerator {
public: // constructors
    MotifStateEnumerator(
            motif_search::SolutionToplogy);

    ~MotifStateEnumerator() {}

public: // main interface
    void
    start(
            structure::BasepairStateOP);

    bool
    finished();

    void
    next();

    motif::MotifStateOP
    top_state();

    motif::MotifStateOPs const &
    all_states();


public: // setters
    inline
    void
    set_size_limit(
            int size_limit) {
        size_limit_ = size_limit;
    }

private:
    void
    _update_current_states();

    void
    _iterate();

    bool
    _within_size_limit();

private:
    Ints indices_;
    Ints maxes_;
    std::vector<motif::MotifStateOPs> motif_states_;
    motif::MotifStateOPs current_;
    structure::BasepairStateOP start_bp_;
    motif::MotifStateAligner aligner_;
    int update_;
    bool updated_;
    int end_;
    int size_limit_;

};

}
}

#endif //RNAMAKE_NEW_MOTIF_STATE_ENUMERATOR_H
