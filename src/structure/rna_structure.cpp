//
//  rna_structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <memory>

//RNAMake Headers
#include "structure/close_chain.h"
#include "structure/rna_structure.h"

namespace structure {

// get basepair functions //////////////////////////////////////////////////////////////////////////

BasepairOPs
structure::RNAStructure::get_basepair(util::Uuid const & bp_uuid) {
    BasepairOPs bps;
    for (auto const & bp : basepairs_) {
        if (bp->uuid() == bp_uuid) { bps.push_back(bp); }
        if (bp->res1()->uuid() == bp_uuid || bp->res2()->uuid() == bp_uuid) { bps.push_back(bp); }

    }
    return bps;
}

BasepairOPs
structure::RNAStructure::get_basepair(
        ResidueOP const & res1,
        ResidueOP const & res2) {
    BasepairOPs bps;
    for (auto & bp : basepairs_) {
        if (bp->res1()->uuid() == res1->uuid() && bp->res2()->uuid() == res2->uuid()) { bps.push_back(bp); }
        if (bp->res1()->uuid() == res2->uuid() && bp->res2()->uuid() == res1->uuid()) { bps.push_back(bp); }
    }
    return bps;
}

BasepairOPs
structure::RNAStructure::get_basepair(
        util::Uuid const & uuid1,
        util::Uuid const & uuid2) {

    BasepairOPs bps;
    for (auto const & bp : basepairs_) {
        if (bp->res1()->uuid() == uuid1 && bp->res2()->uuid() == uuid2) { bps.push_back(bp); }
        if (bp->res1()->uuid() == uuid2 && bp->res2()->uuid() == uuid1) { bps.push_back(bp); }
    }
    return bps;
}

BasepairOPs
structure::RNAStructure::get_basepair(
        String const & name) {
    Strings name_spl = base::split_str_by_delimiter(name, "-");
    String alt_name = name_spl[1] + "-" + name_spl[0];
    for (auto const & bp : basepairs_) {
        if (name.compare(bp->name()) == 0 || alt_name.compare(bp->name()) == 0) {
            return BasepairOPs{bp};
        }
    }

    throw std::runtime_error("could not find basepair with name " + name);

}


// get beads functions /////////////////////////////////////////////////////////////////////////////


util::Beads const &
structure::RNAStructure::get_beads(
        BasepairOPs const & excluded_ends) {

    ResidueOPs excluded_res;
    for (auto const & end : excluded_ends) {
        excluded_res.push_back(end->res1());
        excluded_res.push_back(end->res2());
    }
    beads_ = structure_->get_beads(excluded_res);
    return beads_;
}

util::Beads const &
structure::RNAStructure::get_beads(
        BasepairOP const & excluded_end) {

    ResidueOPs excluded_res;
    excluded_res.push_back(excluded_end->res1());
    excluded_res.push_back(excluded_end->res2());
    beads_ = structure_->get_beads(excluded_res);
    return beads_;
}


// get end infomation functions ////////////////////////////////////////////////////////////////////


int
structure::RNAStructure::get_end_index(BasepairOP const & end) {
    if (std::find(ends_.begin(), ends_.end(), end) == ends_.end()) {
        throw structure::RNAStructureException("cannot find end: " + end->get_name_str() + " in ends ");
    }

    int pos = (int) (std::find(ends_.begin(), ends_.end(), end) - ends_.begin());
    return pos;
}

int
structure::RNAStructure::get_end_index(String const & end_name) {
    for (int i = 0; i < ends_.size(); i++) {
        if (ends_[i]->get_name_str() == end_name) { return i; }
    }

    throw structure::RNAStructureException("could not find a end basepair with name: " + end_name);

    return -1;
}


// output functions ////////////////////////////////////////////////////////////////////////////////


String const
structure::RNAStructure::to_pdb_str(
        int renumber,
        int close_chains) {
    return structure_->to_pdb_str(renumber);
}

void
structure::RNAStructure::to_pdb(
        String const fname,
        int renumber,
        int close_chains,
        int conect_statements) {

    if (close_chains) {
        for (auto & c : structure_->chains()) {
            close_chain(c);
        }
    }

    return structure_->to_pdb(fname, renumber, conect_statements);
}


// non member functions ////////////////////////////////////////////////////////////////////////////


std::shared_ptr<BasepairOPs>
end_from_basepairs(
        StructureOP const & s,
        BasepairOPs const & bps) {

    auto chain_ends = ResidueOPs();
    for (auto const & c : s->chains()) {
        chain_ends.push_back(c->first());
        if (c->length() > 1) {
            chain_ends.push_back(c->last());
        }
    }

    auto res_map = std::map<util::Uuid, ResidueOP, util::UuidCompare>();
    for (auto const & r : chain_ends) { res_map[r->get_uuid()] = r; }

    auto ends = std::make_shared<BasepairOPs>();
    for (auto const & bp : bps) {
        if (bp->get_bp_type() != "cW-W") { continue; }
        if (res_map.find(bp->res1()->uuid()) != res_map.end() &&
            res_map.find(bp->res2()->uuid()) != res_map.end()) {
            ends->push_back(bp);
        }
    }

    return ends;

}

std::shared_ptr<BasepairOPs>
subselect_basepairs_with_res(
        ResidueOPs const & res,
        BasepairOPs const & all_bps) {

    auto res_map = std::map<util::Uuid, ResidueOP, util::UuidCompare>();
    auto bps = std::make_shared<BasepairOPs>();

    for (auto const & r : res) { res_map[r->uuid()] = r; }

    for (auto const & bp : all_bps) {
        if (res_map.find(bp->res1()->uuid()) != res_map.end() &&
            res_map.find(bp->res2()->uuid()) != res_map.end()) {
            bps->push_back(bp);
        }
    }

    return bps;
}

}












