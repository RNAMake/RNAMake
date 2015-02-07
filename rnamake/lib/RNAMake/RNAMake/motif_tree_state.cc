//
//  motif_tree_state.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state.h"
#include "types.h"
#include "FileIO.h"

NameElements
parse_db_name(
    String const & s) {
    Strings spl = split_str_by_delimiter(s, "-");
    return NameElements(spl[0], std::stoi(spl[1]), std::stoi(spl[2]),
                        std::stoi(spl[3]), std::stoi(spl[4]), std::stoi(spl[5]), std::stoi(spl[6]));
}

MotifTreeStateOP
motif_to_state(
    MotifOP const & m,
    int end_index,
    int end_flip) {
    
    MotifTree mt;
    mt.add_motif(m);
    MotifOP m_copy = mt.nodes()[1]->motif();
    std::stringstream ss;
    ss << m_copy->name() << "-" << end_index << "-" << end_flip;
    BasepairStateOPs ends ( m_copy->ends().size());
    int i = -1;
    for (auto const & end : m_copy->ends()) {
        i++;
        if(i == end_index) { continue; }
        ends[i] = BasepairStateOP(new BasepairState(end->state()));
    }
    Points beads;
    for(auto const & b : m_copy->beads()) {
        if (b.btype() != PHOS) { beads.push_back(b.center()); }
    }
    String build_str = m_copy->to_str();
    String name = ss.str();
    float score = 0.0;
    MotifTreeStateOP mts ( new MotifTreeState(name, end_index, (int)m_copy->residues().size(), score, beads, ends, end_flip,  build_str));
    return mts;
}