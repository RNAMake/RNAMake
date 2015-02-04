//
//  motif.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif.h"
#include "xyzMatrix.h"
#include "transform.h"

Motif::Motif(
    String const & s,
    ResidueTypeSet const & rts):
    beads_(Beads()),
    score_(0),
    basepairs_(BasepairOPs()),
    ends_(BasepairOPs()),
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
        ResidueOP res1 = structure_.get_residue(res1_num, res1_id, "");
        ResidueOP res2 = structure_.get_residue(res2_num, res2_id, "");
        BasepairState bpstate = str_to_basepairstate(bp_spl[1]);
        BasepairOP bp ( new Basepair(res1, res2, bpstate.r(), bp_spl[1] ));
        bp->flip( std::stoi(bp_spl[4]));
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
    cmotif.basepairs_ = BasepairOPs();
    cmotif.ends_ = BasepairOPs();
    int i = 0;
    for (auto const & b : beads_) {
        cmotif.beads_[i] = b.copy();
        i++;
    }
    i = 0;
    //ResidueOP res1, res2;
    /*std::cout << cmotif.residues().size() << std::endl;
    for (auto const & r : residues()) {
        std::cout << r->num() << " " << r->uuid().s_uuid() << std::endl;
    }*/
    
    //std::cout << basepairs_[0].res1().num() << std::endl;
    //exit(0);
    for (auto const & bp : basepairs_) {
        ResidueOP res1 = cmotif.get_residue(bp->res1()->uuid());
        ResidueOP res2 = cmotif.get_residue(bp->res2()->uuid());
        BasepairOP new_bp ( new Basepair ( res1, res2, bp->r(), bp->bp_type() )) ;
        new_bp->flip(bp->flipped());
        new_bp->uuid(bp->uuid());
        cmotif.basepairs_.push_back(new_bp);
    }
    i = 0;
    for (auto & end: ends_) {
        BasepairOPs bps = cmotif.get_basepair(end->uuid());
        cmotif.ends_.push_back(bps[0]);
    }
    
    cmotif._cache_basepair_frames();
    return cmotif;
}

String const
Motif::to_str() {
    std::stringstream ss;
    ss << mdir_ << "&" << name_ << "&" << score_ << "&" << mtype_ << "&" << structure_.to_str() << "&";
    for ( auto const & bp : basepairs_ ) {
        ss << bp->to_str() << "@";
    }
    ss << "&";
    for ( auto const & end : ends_) {
        int pos = 0;
        int i = 0;
        for (auto const & bp : basepairs_) {
            if( end == bp) { pos = i; break; }
            i ++;
        }
        //(int)(std::find(basepairs_.begin(), basepairs_.end(), end) - basepairs_.begin());
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


BasepairOPs
Motif::get_basepair(Uuid const & bp_uuid) {
    BasepairOPs bps;
    for( auto const & bp : basepairs_) {
        if(bp->uuid() == bp_uuid) { bps.push_back(bp); }
    }
    return bps;
}

BasepairOPs
Motif::get_basepair(ResidueOP res1,
                    ResidueOP res2) {
    BasepairOPs bps;
    for(auto & bp : basepairs_) {
        if(bp->res1()->uuid() == res1->uuid() && bp->res2()->uuid() == res2->uuid()) { bps.push_back(bp); }
        if(bp->res1()->uuid() == res2->uuid() && bp->res2()->uuid() == res1->uuid()) { bps.push_back(bp); }
    }
    return bps;
}

BasepairOPs
Motif::get_basepair(Uuid const & uuid1,
                    Uuid const & uuid2) {
    BasepairOPs bps;
    for( auto const & bp : basepairs_) {
        if(bp->res1()->uuid() == uuid1 && bp->res2()->uuid() == uuid2) { bps.push_back(bp); }
        if(bp->res1()->uuid() == uuid2 && bp->res2()->uuid() == uuid1) { bps.push_back(bp); }
    }
    return bps;
}


void
align_motif(BasepairOP ref_bp,
            BasepairOP motif_end,
            Motif & motif) {
    
    Matrix ref_T;
    transpose(ref_bp->r(), ref_T);
    Matrix r;
    dot(ref_T, motif_end->r(), r);
    Point trans = -motif_end->d();
    Transform t(r, trans);
    motif.transform(t);
    Point bp_pos_diff = ref_bp->d() - motif_end->d();
    std::cout << bp_pos_diff << std::endl;
    motif.move(bp_pos_diff);
    bp_pos_diff = ref_bp->d() - motif_end->d();
    std::cout << bp_pos_diff << std::endl;

}


