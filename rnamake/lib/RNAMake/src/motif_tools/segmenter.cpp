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
    RNAStructureOP const & m,
    ResidueOPs const & res) {
    
    pairs_ = PairOPs();
    end_pairs_ = PairOPs();
    
    ChainOP sc, sc1, sc2;
    float dist;
    
    for(auto const & c : m->chains()) {
        int i = -1;
        for(auto const & res1 : res) {
            i++;
            int j = -1;
            for(auto const & res2 : res) {
                j++;
                if(i >= j) { continue; }
                try {
                    sc = c->subchain(res1, res2);
                } catch(ChainException) {continue; }
                dist = (int)sc->length();
                pairs_.push_back(std::make_shared<Pair>(res1, res2, dist));
            }
            try {
                sc1 = c->subchain(res1, c->last());
            } catch(ChainException) { continue; }
            try {
                sc2 = c->subchain(c->first(), res1);
            } catch(ChainException) { continue; }
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
    RNAStructureOP const & m,
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
    RNAStructureOP const & m,
    ResidueOPs & res,
    BasepairOPs const & bps,
    ResidueOPs const & cutpoints,
    BasepairOPs const & cut_bps) {

    auto other_res = ResidueOPs();
    for(auto const & r : cutpoints) { other_res.push_back(r); }
    for(auto const & r : m->residues()) {
        int found = 0;
        for(auto const & r2 : res) {
            if(r->uuid() == r2->uuid()) {
                found = 1; break;
            }
        }
        if(found == 0) { other_res.push_back(r); }
        
    }
    
    auto res_uuids = std::map<Uuid, int, UuidCompare>();
    
    for(auto const & r : other_res) {
        res_uuids[r->uuid()] = 1;
    }
    
    auto removed = mf_.motif_from_res(res, bps);
    removed->name(m->name() + ".removed");
    mf_.standardize_rna_structure_ends(removed);

    auto other_bps = BasepairOPs();
    
    for(auto const & bp : m->basepairs()) {
        if(res_uuids.find(bp->res1()->uuid()) != res_uuids.end() &&
           res_uuids.find(bp->res2()->uuid()) != res_uuids.end()) {
            other_bps.push_back(bp);
        }
    }
    
    auto remaining = mf_.motif_from_res(other_res, other_bps);
    remaining->name(m->name() + ".remaining");
    mf_.standardize_rna_structure_ends(remaining);
    
    for(auto & end : remaining->ends()) {
        int flip_res = 0;
        for (auto const & c : remaining->chains()) {
            if(c->first() == end->res2()) {
                flip_res = 1;
                break;
            }
        }
        
        if(!flip_res) { continue; }
        
        auto temp = end->res1();
        end->res1(end->res2());
        end->res2(temp);
    }
    
    
    return std::make_shared<Segments>(removed, remaining);
    
}


SegmentsOP
Segmenter::apply(
    RNAStructureOP const & m,
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
        return _get_segments(m, sub_res_array, basepairs, res, bps);
    }
    
    throw SegmenterException("failed to segment rna structure");
}












