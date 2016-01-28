//
//  segmenter.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__segmenter__
#define __RNAMake__segmenter__

#include <stdio.h>
#include "structure/residue.h"
#include "motif/pose.h"

struct Pair {
    Pair(
        ResidueOP const & nres1,
        ResidueOP const & nres2,
        int ndist):
    res1(nres1),
    res2(nres2),
    dist(ndist)
    {}
    
    int
    inline
    contains(
        ResidueOP const & res) {
        
        if(res == res1 || res == res2) { return 1; }
        else                           { return 0; }
    }
    
    ResidueOP res1, res2;
    int dist;
};

struct Segments {
    Segments(
        PoseOP const & nremaining,
        PoseOP const & nremoved):
    removed(nremoved),
    remaining(nremaining)
    {}
    
    PoseOP removed, remaining;
};

typedef std::shared_ptr<Pair> PairOP;
typedef std::vector<PairOP>   PairOPs;
typedef std::shared_ptr<Segments> SegmentsOP;

struct PairSearchNode {
    PairSearchNode(
        PairOPs const & npairs):
    pairs(npairs) {}
    
    ~PairSearchNode() {}
    
    inline
    int
    contains(
        ResidueOP const & res) {
        
        for(auto const & p : pairs) {
            if(p->contains(res)) { return 1; }
        }
        return 0;
    }
    
    PairOPs pairs;
    int score;
};


class Segmenter {
public:
    Segmenter() {}
    
    ~Segmenter() {}
    
public:
    SegmentsOP
    apply(
        PoseOP const &,
        BasepairOPs const &);
    
private:
    
    void
    _get_pairs(
        PoseOP const &,
        ResidueOPs const &);
    
private:
    
    PairOPs pairs_, end_pairs_;
    
};

#endif /* defined(__RNAMake__segmenter__) */
