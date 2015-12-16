//
//  rna_structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <memory>

//RNAMake Headers
#include "structure/rna_structure.h"

BasepairOPs
RNAStructure::get_basepair(Uuid const & bp_uuid) {
    BasepairOPs bps;
    for( auto const & bp : basepairs_) {
        if(bp->uuid() == bp_uuid) { bps.push_back(bp); }
        if(bp->res1()->uuid() == bp_uuid || bp->res2()->uuid() == bp_uuid ) { bps.push_back(bp); }
        
    }
    return bps;
}

BasepairOPs
RNAStructure::get_basepair(
    ResidueOP const & res1,
    ResidueOP const & res2) {
    BasepairOPs bps;
    for(auto & bp : basepairs_) {
        if(bp->res1()->uuid() == res1->uuid() && bp->res2()->uuid() == res2->uuid()) { bps.push_back(bp); }
        if(bp->res1()->uuid() == res2->uuid() && bp->res2()->uuid() == res1->uuid()) { bps.push_back(bp); }
    }
    return bps;
}

BasepairOPs
RNAStructure::get_basepair(
    Uuid const & uuid1,
    Uuid const & uuid2) {
    
    BasepairOPs bps;
    for( auto const & bp : basepairs_) {
        if(bp->res1()->uuid() == uuid1 && bp->res2()->uuid() == uuid2) { bps.push_back(bp); }
        if(bp->res1()->uuid() == uuid2 && bp->res2()->uuid() == uuid1) { bps.push_back(bp); }
    }
    return bps;
}

Beads const &
RNAStructure::get_beads(
    BasepairOPs const & excluded_ends) {
    
    ResidueOPs excluded_res;
    for (auto const & end : excluded_ends) {
        excluded_res.push_back(end->res1());
        excluded_res.push_back(end->res2());
    }
    beads_ = structure_->get_beads(excluded_res);
    return beads_;
}

Beads const &
RNAStructure::get_beads(
    BasepairOP const & excluded_end) {
    
    ResidueOPs excluded_res;
    excluded_res.push_back(excluded_end->res1());
    excluded_res.push_back(excluded_end->res2());
    beads_ = structure_->get_beads(excluded_res);
    return beads_;
}

BasepairOPs
RNAStructure::get_basepair(
    String const & name) {
    Strings name_spl = split_str_by_delimiter(name, "-");
    String alt_name = name_spl[1] + "-" + name_spl[0];
    for(auto const & bp : basepairs_) {
        if(name.compare(bp->name()) == 0 || alt_name.compare(bp->name()) == 0) {
            return BasepairOPs{ bp };
        }
    }
    
    throw "could not find basepair with name " + name;
    
}

int
RNAStructure::end_index(BasepairOP const & end) {
    int pos = (int)(std::find(ends_.begin(), ends_.end(), end) - ends_.begin());
    return pos;
}

String const
RNAStructure::to_pdb_str(
    int renumber) {
    return structure_->to_pdb_str(renumber);
}

void
RNAStructure::to_pdb(
    String const fname,
    int renumber) {
    return structure_->to_pdb(fname, renumber);
}


std::unique_ptr<BasepairOPs>
end_from_basepairs(
    StructureOP const & s,
    BasepairOPs const & bps) {
    
    auto chain_ends = ResidueOPs();
    for(auto const & c : s->chains()) {
        chain_ends.push_back(c->first());
        if(c->length() > 1) {
            chain_ends.push_back(c->last());
        }
    }
    
    auto res_map = std::map<Uuid, ResidueOP, UuidCompare>();
    for(auto const & r : chain_ends ) { res_map[r->uuid()] = r;}

    auto ends = std::make_unique<BasepairOPs>();
    for(auto const & bp : bps) {
        if(bp->bp_type() != "cW-W") { continue; }
        if(res_map.find(bp->res1()->uuid()) != res_map.end() &&
           res_map.find(bp->res2()->uuid()) != res_map.end() ) {
            ends->push_back(bp);
        }
    }
    
    return ends;
    
}

std::unique_ptr<BasepairOPs>
subselect_basepairs_with_res(
    ResidueOPs const & res,
    BasepairOPs const & all_bps) {
    
    auto res_map = std::map<Uuid, ResidueOP, UuidCompare>();
    auto bps = std::make_unique<BasepairOPs>();
    
    for(auto const & r : res ) { res_map[r->uuid()] = r;}

    for(auto const & bp : all_bps) {
        if(res_map.find(bp->res1()->uuid()) != res_map.end() &&
           res_map.find(bp->res2()->uuid()) != res_map.end() )  {
            bps->push_back(bp);
        }
    }

    return bps;
}












