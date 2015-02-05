//
//  motif_tree_state_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_library.h"
#include "motif_type.h"
#include "settings.h"
#include "types.h"
#include "xyzVector.h"
#include "basepair_state.h"

MotifTreeStateLibrary::MotifTreeStateLibrary(MotifType const & mtype) {
    mtype_ = mtype;
    String path = resources_path() + "/precomputed/motif_tree_states/" + type_to_str(mtype_) + ".new.me";
    _load_states_from_file(path);
    
}

void
MotifTreeStateLibrary::_load_states_from_file(String const & path) {
    String line;
    std::ifstream in;
    motif_tree_states_ = MotifTreeStates();
    in.open(path.c_str());
    String name, build_string;
    float score;
    int size, flip;
    Points beads;
    while ( in.good() ) {
        getline(in, line);
        if(line.length() < 10) { continue; }
        Strings spl = split_str_by_delimiter(line, "|");
        name = spl[0];
        score = std::stof(spl[1]);
        size = std::stoi(spl[2]);
        flip = std::stoi(spl[3]);
        build_string = spl[4];
        beads = vectors_from_str(spl[5]);
        BasepairStateOPs end_states;
        for (int i = 6; i < spl.size(); i++) {
            if(spl[i].length() < 5) { end_states.push_back(NULL); }
            else {
                BasepairStateOP bpstate ( new BasepairState ( str_to_basepairstate(spl[i])));
                end_states.push_back(bpstate);
            }
        }
        NameElements name_elements = parse_db_name(name);
        MotifTreeState motif_tree_state ( name, name_elements.start_index, size, score, beads, end_states, flip,
                                         build_string);
        motif_tree_states_.push_back(motif_tree_state);
    }
}
