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
#include <queue>
#include <map>
#include "structure/residue.h"
#include "motif/motif.h"

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
    MotifOP removed, remaining;
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

typedef std::vector<PairSearchNode> PairSearchNodes;

struct PairSearchNodeCompare {
    bool
    operator () (
        PairSearchNode const & node1,
        PairSearchNode const & node2) {
        
        if (node1.score > node2.score) { return true;  }
        else     					   { return false; }
    }
};

typedef std::priority_queue<PairSearchNode, PairSearchNodes, PairSearchNodeCompare> PairSearchNodePriortyQueue;

class PairSearch {
public:
    PairSearch():
    queue_(PairSearchNodePriortyQueue()),
    solutions_(PairSearchNodes()),
    values_(std::map<ResidueOP, float>())
    {}
    
    ~PairSearch() {}
    
    
public:
    
    PairSearchNodes const &
    search(
        ResidueOPs const & res,
        PairOPs const & pairs,
        PairOPs const & end_pairs) {
        _get_default_values(res, pairs, end_pairs);
    
        
        return solutions_;
    
    }
    
private:
    void
    _get_default_values(
        ResidueOPs const & res,
        PairOPs const & pairs,
        PairOPs const & end_pairs) {
        
        auto all_pairs = PairOPs();
        for(auto const & p : pairs) { all_pairs.push_back(p); }
        for(auto const & p : end_pairs) { all_pairs.push_back(p); }
        
        for(auto const & r : res) {
            auto total = 0.0, count = 0.0;
            for(auto const & p : all_pairs) {
                if(p->contains(r)) {
                    total += p->dist;
                    count += 1;
                }
            }
            values_[r] = (total / count) / 2;
        }
    }
    
  
    
private:
    PairSearchNodePriortyQueue queue_;
    PairSearchNodes solutions_;
    std::map<ResidueOP, float> values_;
};


class Segmenter {
public:
    Segmenter() {}
    
    ~Segmenter() {}
    
public:
    SegmentsOP
    apply(
        MotifOP const &,
        BasepairOPs const &);
    
private:
    
    void
    _get_pairs(
        MotifOP const &,
        ResidueOPs const &);
    
private:
    
    PairOPs pairs_, end_pairs_;
    
};


#endif /* defined(__RNAMake__segmenter__) */
