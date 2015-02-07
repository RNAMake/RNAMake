//
//  motif_tree_state.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state__
#define __RNAMake__motif_tree_state__

#include <stdio.h>
#include "types.h"
#include "xyzVector.h"
#include "basepair_state.h"
#include "motif_tree.h"

struct NameElements {
public:
    NameElements(
        String const & nmotif_name,
        int const & nhelix_direction,
        int const & nstart_helix_count,
        int const & nstart_index,
        int const & nend_helix_count,
        int const & nend_index,
        int const & nflip_direction):
        motif_name ( nmotif_name),
        helix_direction( nhelix_direction),
        start_helix_count( nstart_helix_count),
        start_index( nstart_index),
        end_helix_count( nend_helix_count),
        end_index( nend_index),
        flip_direction( nflip_direction)
    {}
    
    ~NameElements() {}

public:
    String motif_name;
    int helix_direction, start_helix_count, start_index, end_helix_count, end_index, flip_direction;
    
    
};

class MotifTreeState {
public:
    MotifTreeState(
        String const & name,
        int const & start_index,
        int const & size,
        float const & score,
        Points const & beads,
        BasepairStateOPs const & ends,
        int const & flip,
        String const & build_string):
        name_( name),
        start_index_(start_index),
        size_(size),
        score_(score),
        beads_(beads),
        ends_(ends),
        flip_(flip),
        build_string_(build_string)
    {}
    
    ~MotifTreeState() {}
    
public: //getters
    
    inline
    String const &
    name() const { return name_; }
    
    inline
    String const &
    build_string() const { return build_string_; }

    inline
    int const &
    start_index() const { return start_index_; }

    inline
    int const &
    size() const { return size_; }
    
    inline
    int const &
    flip() const { return flip_; }
   
    inline
    float const &
    score() const { return score_; }
    
    inline
    Points const &
    beads() const { return beads_; }
    
    inline
    BasepairStateOPs const &
    end_states() const { return ends_; }
    
private:
    String name_, build_string_;
    int start_index_, size_, flip_;
    float score_;
    Points beads_;
    BasepairStateOPs ends_;
};


NameElements
parse_db_name(
    String const &);


typedef std::vector<MotifTreeState> MotifTreeStates;
typedef std::shared_ptr<MotifTreeState> MotifTreeStateOP;
typedef std::vector<MotifTreeStateOP> MotifTreeStateOPs;

MotifTreeStateOP
motif_to_state(MotifOP const &,
               int end_index = 0,
               int end_flip = 0);


#endif /* defined(__RNAMake__motif_tree_state__) */
