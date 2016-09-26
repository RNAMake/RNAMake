//
//  segmenter.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <algorithm>

#include "structure/chain.h"
#include "motif_tools/segmenter.h"

void
Segmenter::_get_pairs(
    MotifOP const & m,
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
            }
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


ChainOP
Segmenter::_get_subchain(
    MotifOP const & m,
    PairOP const & pair) {
    
    for(auto const & c : m->chains()) {
        try {
            auto sc = c->subchain(pair->res1, pair->res2);
            return sc;
        }
        catch(ChainException) { continue; }
    }
    
    throw SegmenterException("cannot get sub chain");
    
}


SegmentsOP
Segmenter::_get_segments(
    MotifOP const & m,
    ResidueOPs & res,
    BasepairOPs const & bps,
    ResidueOPs const & cutpoints) {

    auto removed = mf_.motif_from_res(res, bps);
    removed->name(m->name() + ".removed");

    auto other_res = ResidueOPs();
    for(auto const & r : cutpoints) { other_res.push_back(r); }
    auto other_bps = BasepairOPs();
    
    for(auto const & r : m->residues()) {
        if(std::find(res.begin(), res.end(), r) == res.end()) {
            other_res.push_back(r);
        }
    }
    
    for(auto const & bp : m->basepairs()) {
        if(std::find(bps.begin(), bps.end(), bp) == bps.end()) {
            other_bps.push_back(bp);
        }
    }

    auto remaining = MotifOP();
    try {
        remaining = mf_.motif_from_res(other_res, other_bps);
        remaining->name(m->name() + ".remaining");
    }
    catch(MotifFactoryException const & e) {
        if(res.size() > 6) { throw e; }
    }
    
    return std::make_shared<Segments>(Segments{removed, remaining});
    
}


SegmentsOP
Segmenter::apply(
    MotifOP const & m,
    BasepairOPs const & bps) {
    
    ResidueOPs res;
    for(auto const & bp : bps) {
        for(auto const & r : bp->residues()) {
            res.push_back(r);
        }
    }
    
    _get_pairs(m, res);
    auto pair_search = PairSearch();
    auto sols = pair_search.search(res, pairs_, end_pairs_);

    for(auto const & s : sols) {
        auto subchains = ChainOPs();
        auto sub_res = std::map<Uuid, int, UuidCompare>();
        auto sub_res_array = ResidueOPs();
        
        for(auto const & r: res) {
            sub_res[r->uuid()] = 1;
            sub_res_array.push_back(r);
        }
        
        for(auto const & p : s.pairs) {
            auto sc = _get_subchain(m, p);
            subchains.push_back(sc);
            for(auto const & r : sc->residues()) {
                sub_res[r->uuid()] = 1;
                sub_res_array.push_back(r);
            }
        }
        
        auto basepairs = BasepairOPs();
        int missed_bps = 0;
        
        for(auto const & bp : m->basepairs()) {
            if(bp->bp_type() != "cW-W") { continue; }
            if(sub_res.find(bp->res1()->uuid()) != sub_res.end() &&
               sub_res.find(bp->res2()->uuid()) != sub_res.end()) {
                basepairs.push_back(bp);
            }
            else if(sub_res.find(bp->res1()->uuid()) != sub_res.end()) { missed_bps += 1; }
            else if(sub_res.find(bp->res2()->uuid()) != sub_res.end()) { missed_bps += 1; }
        }
        
        if(missed_bps > 2) { continue; }
        return _get_segments(m, sub_res_array, basepairs, res);
    }
    
    return nullptr;
}












