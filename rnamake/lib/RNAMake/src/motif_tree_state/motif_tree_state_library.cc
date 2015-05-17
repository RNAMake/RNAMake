//
//  motif_tree_state_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/types.h"
#include "util/settings.h"
#include "util/file_io.h"
#include "math/xyz_vector.h"
#include "structure/basepair_state.h"
#include "motif/motif_type.h"
#include "motif_tree_state/motif_tree_state_library.h"

MotifTreeStateLibrary::MotifTreeStateLibrary(
    MotifType const & mtype,
    int new_s) {
    mtype_ = mtype;
    new_s_ = new_s;
    String path = resources_path() + "/precomputed/motif_tree_states/" + type_to_str(mtype_) + ".new.me";
    _load_states_from_file(path, 0);
    null_ = nullptr;
    
}

MotifTreeStateLibrary::MotifTreeStateLibrary(MotifType const & mtype, int size_limit, int new_s) {
    mtype_ = mtype;
    String path = resources_path() + "/precomputed/motif_tree_states/" + type_to_str(mtype_) + ".new.me";
    new_s_ = new_s;
    _load_states_from_file(path, size_limit);
    null_ = nullptr;
    
}

MotifTreeStateLibrary::MotifTreeStateLibrary(String const & libpath, int new_s) {
    mtype_ = UNKNOWN;
    new_s_ = new_s;
    _load_states_from_file(libpath, 0);
    null_ = nullptr;
}

void
MotifTreeStateLibrary::_load_states_from_file(String const & path, int size_limit) {
    motif_tree_states_ = MotifTreeStateOPs();
    motif_tree_state_dict_ = std::map<String, MotifTreeStateOP>();

    Strings lines = get_lines_from_file(path);
    
    String name, build_string;
    float score;
    int size, flip;
    Points beads;
    for(auto const & line : lines ) {
        if(line.length() < 10) { continue; }
        Strings spl = split_str_by_delimiter(line, "|");
        name = spl[0];
        score = std::stof(spl[1]);
        size = std::stoi(spl[2]);
        if(size_limit > size) { continue; }
        flip = std::stoi(spl[3]);
        build_string = spl[4];
        beads = vectors_from_str(spl[5]);
        BasepairStateOPs end_states;
        for (int i = 6; i < spl.size(); i++) {
            if(spl[i].length() < 5) { end_states.push_back(nullptr); }
            else {
                BasepairStateOP bpstate ( new BasepairState ( str_to_basepairstate(spl[i])));
                end_states.push_back(bpstate);
            }
        }
        if(new_s_ == 0) {
            NameElements name_elements = parse_db_name(name);
            MotifTreeStateOP motif_tree_state (new MotifTreeState ( name, name_elements.start_index, size, score, beads, end_states, flip, build_string));
            motif_tree_states_.push_back(motif_tree_state);
            motif_tree_state_dict_[name] = motif_tree_state;
        }
        else {
            Strings spl = split_str_by_delimiter(name, "-");
            int start_index = std::stoi(spl[1]);
            MotifTreeStateOP motif_tree_state (new MotifTreeState ( name, start_index, size, score, beads, end_states, flip, build_string));
            motif_tree_states_.push_back(motif_tree_state);
            motif_tree_state_dict_[name] = motif_tree_state;
        }
    }
}

MotifTreeStateOP const &
MotifTreeStateLibrary::get_state(String const & name) {
    for(auto const & mts : motif_tree_states_) {
        if(mts->name().compare(name) == 0) { return mts; }
    }
    throw "could not find mts in library";
}

MotifTreeStateOP const &
MotifTreeStateLibrary::get_state_no_error(String const & name) {
    for(auto const & mts : motif_tree_states_) {
        if(mts->name().compare(name) == 0) { return mts; }
    }
    return null_;
}



