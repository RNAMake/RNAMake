//
//  segmenter.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "structure/chain.h"
#include "motif_tools/segmenter.h"

/*
void
Segmenter::_get_pairs(
    PoseOP const & m,
    ResidueOPs const & res) {
    
    pairs_ = PairOPs();
    end_pairs_ = PairOPs();
    
    ChainOP sc, sc1, sc2;
    int i = 0, j = 0;
    float dist;
    
    for(auto const & c : m->chains()) {
        i = -1;
        for(auto const & res1 : res) {
            i++;
            j = -1;
            for(auto const & res2 : res) {
                j++;
                if(i >= j) { continue; }
                sc = c->subchain(res1, res2);
                if(sc == nullptr) { continue; }
                dist = (int)sc->length();
                pairs_.push_back(std::make_shared<Pair>(res1, res2, dist));
                sc1 = c->subchain(res1, c->last());
                sc2 = c->subchain(c->first(), res1);
                if(sc1 == nullptr) { continue; }
                if(res1 != c->last() &&
                   std::find(res.begin(), res.end(), c->last()) == res.end()) {
                    auto pair = std::make_shared<Pair>(res1, c->last(), sc1->length());
                    end_pairs_.push_back(pair);
                }
                if(res1 != c->first() &&
                   std::find(res.begin(), res.end(), c->first()) == res.end()) {
                    auto pair = std::make_shared<Pair>(res1, c->first(), sc2->length());
                    end_pairs_.push_back(pair);
                }
            }
        }
    }
}
*/

SegmentsOP
Segmenter::apply(
    PoseOP const & m,
    BasepairOPs const & bps) {
    
    /*ResidueOPs res;
    for(auto const & bp : bps) {
        for(auto const & r : bp->residues()) {
            res.push_back(r);
        }
    }
    
    _get_pairs(m, res);
    */
    return nullptr;
}
