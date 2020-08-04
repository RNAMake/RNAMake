//
//  motif_state_node.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef motif_state_node_hpp
#define motif_state_node_hpp

#include <stdio.h>
#include "motif/motif_state.h"

namespace motif_data_structure {

struct MSNodeData {
public:
    inline
    MSNodeData(
            motif::MotifStateOP const & nref_state) :
            ref_state(nref_state),
            cur_state(std::make_shared<motif::MotifState>(*nref_state)) {}

    inline
    MSNodeData(
            MSNodeData const & ndata) :
            ref_state(std::make_shared<motif::MotifState>(*ndata.ref_state)),
            cur_state(std::make_shared<motif::MotifState>(*ndata.cur_state)) {}

public: //wrappers for current state

    inline
    structure::BasepairStateOP
    get_end_state(String const & name) { return cur_state->get_end_state(name); }

    inline
    structure::BasepairStateOP
    get_end_state(int i) { return cur_state->end_states()[i]; }

    inline
    int
    get_end_index(String const & name) { return cur_state->get_end_index(name); }

    inline
    String const &
    name() { return cur_state->name(); }

    inline
    int
    block_end_add() { return cur_state->block_end_add(); }

    inline
    String const &
    end_name(int i) { return cur_state->end_names()[i]; }

    inline
    util::Uuid const &
    uuid() { return cur_state->uuid(); }

public: //wrappers to set some values

    inline
    void
    uuid(util::Uuid const & uuid) {
        cur_state->uuid(uuid);
        ref_state->uuid(uuid);
    }


public:
    motif::MotifStateOP ref_state, cur_state;

};

typedef std::shared_ptr<MSNodeData> MSNodeDataOP;

}

#endif /* motif_state_node_hpp */
