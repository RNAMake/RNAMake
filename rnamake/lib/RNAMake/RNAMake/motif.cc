//
//  motif.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif.h"

Motif::Motif(
    String const & s,
    ResidueTypeSet const & rts):
    beads_(Beads()),
    score_(0),
    basepairs_(Basepairs()),
    ends_(Basepairs()),
    mdir_(String()),
    name_(String()),
    cached_rotations_(Matrices())
{
    
    Strings spl = split_str_by_delimiter(s, "&");
    mdir_ = spl[0];
    name_ = spl[1];
    score_ = std::stof(spl[2]);
    mtype_ = static_cast<MotifType>(std::stoi(spl[3]));
    structure_ = str_to_structure(spl[4], rts);
    Strings basepair_str = split_str_by_delimiter(spl[5], "@");
    for (auto const & bp_str : basepair_str) {
        Strings bp_spl = split_str_by_delimiter(bp_str, ",");
        Strings res_spl = split_str_by_delimiter(bp_spl[0], "-");
        String res1_id = res_spl[0].substr(0,1);
        String res2_id = res_spl[1].substr(0,1);
        int res1_num = std::stoi(res_spl[0].substr(1));
        int res2_num = std::stoi(res_spl[1].substr(1));
        Residue res1 = structure_.get_residue(res1_num, res1_id, "");
        Residue res2 = structure_.get_residue(res2_num, res2_id, "");
        BasepairState bpstate = str_to_basepairstate(bp_spl[1]);
        Basepair bp ( res1, res2, bpstate.r(), bp_spl[1] );
        bp.flip( std::stoi(bp_spl[4]));
        basepairs_.push_back(bp);
    }
    
    Strings end_indexes = split_str_by_delimiter(spl[6], " ");
    for (auto const & index : end_indexes) {
        ends_.push_back( basepairs_ [ std::stoi(index) ]);
    }
    _cache_basepair_frames();
}

Motif
Motif::copy() {
    Motif cmotif;
    cmotif.name_ = name_;
    cmotif.mdir_ = mdir_;
    cmotif.score_ = score_;
    cmotif.mtype_ = mtype_;
    cmotif.structure_ = structure_.copy();
    cmotif.beads_ = Beads(beads_.size());
    cmotif.basepairs_ = Basepairs(basepairs_.size());
    cmotif.ends_ = Basepairs(ends_.size());
    int i = 0;
    for (auto const & b : beads_) {
        cmotif.beads_[i] = b.copy();
        i++;
    }
    i = 0;
    Residue res1, res2;
    for (auto const & bp : basepairs_) {
        res1 = cmotif.get_residue(bp.res1().uuid());
        res2 = cmotif.get_residue(bp.res2().uuid());
        Basepair new_bp ( res1, res2, bp.r(), bp.bp_type() );
        new_bp.flip(bp.flipped());
        new_bp.uuid(bp.uuid());
        cmotif.basepairs_[i] = new_bp;
        i++;
    }
    i = 0;
    for (auto const & end: ends_) {
        Basepairs bps = get_basepair(end.res1(), end.res2());
        cmotif.ends_[i] = bps[0];
        i++;
    }
    
    cmotif._cache_basepair_frames();
    return cmotif;
}

String const
Motif::to_str() {
    std::stringstream ss;
    ss << mdir_ << "&" << name_ << "&" << score_ << "&" << mtype_ << "&" << structure_.to_str() << "&";
    for ( auto const & bp : basepairs_ ) {
        ss << bp.to_str() << "@";
    }
    ss << "&";
    for ( auto const & end : ends_) {
        int pos = (int)(std::find(basepairs_.begin(), basepairs_.end(), end) - basepairs_.begin());
        ss << pos << " ";
    }
    ss << "&";
    return ss.str();
}

String const
Motif::to_pdb_str() {
    return structure_.to_pdb_str();
}

void
Motif::to_pdb(String const fname) {
    return structure_.to_pdb(fname);
}


Basepairs
Motif::get_basepair(Uuid const & bp_uuid) {
    Basepairs bps;
    for( auto const & bp : basepairs_) {
        if(bp.uuid() == bp_uuid) { bps.push_back(bp); }
    }
    return bps;
}

Basepairs
Motif::get_basepair(Residue const & res1,
                    Residue const & res2) {
    Basepairs bps;
    for( auto const & bp : basepairs_) {
        if(bp.res1() == res1 && bp.res2() == res2) { bps.push_back(bp); }
        if(bp.res1() == res2 && bp.res2() == res1) { bps.push_back(bp); }
    }
    return bps;
}

Basepairs
Motif::get_basepair(Uuid const & uuid1,
                    Uuid const & uuid2) {
    Basepairs bps;
    for( auto const & bp : basepairs_) {
        if(bp.res1().uuid() == uuid1 && bp.res2().uuid() == uuid2) { bps.push_back(bp); }
        if(bp.res1().uuid() == uuid2 && bp.res2().uuid() == uuid1) { bps.push_back(bp); }
    }
    return bps;
}






