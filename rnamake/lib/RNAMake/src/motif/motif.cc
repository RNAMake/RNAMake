//
//  motif.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include <assert.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

//RNAMake Headers
#include "util/file_io.h"
#include "util/settings.h"
#include "util/x3dna.h"
#include "structure/residue_type_set_manager.h"
#include "structure/chain.h"
#include "motif/motif.h"

Motif::Motif(
    String const & s,
    ResidueTypeSet const & rts):
    beads_(Beads()),
    score_(0),
    basepairs_(BasepairOPs()),
    ends_(BasepairOPs()),
    path_(String()),
    name_(String())
{
 
    if(s.length() < 10) {
        throw "tried to construct Motif object from string, with a string too short";
    }
    
    Strings spl = split_str_by_delimiter(s, "&");
    path_ = spl[0];
    name_ = spl[1];
    score_ = std::stof(spl[2]);
    block_end_add_ = std::stoi(spl[3]);
    mtype_ = static_cast<MotifType>(std::stoi(spl[4]));
    structure_ = StructureOP( new Structure(str_to_structure(spl[5], rts)));
    Strings basepair_str = split_str_by_delimiter(spl[6], "@");
    for (auto const & bp_str : basepair_str) {
        Strings bp_spl = split_str_by_delimiter(bp_str, ",");
        Strings res_spl = split_str_by_delimiter(bp_spl[0], "-");
        String res1_id = res_spl[0].substr(0,1);
        String res2_id = res_spl[1].substr(0,1);
        int res1_num = std::stoi(res_spl[0].substr(1));
        int res2_num = std::stoi(res_spl[1].substr(1));
        ResidueOP res1 = structure_->get_residue(res1_num, res1_id, "");
        ResidueOP res2 = structure_->get_residue(res2_num, res2_id, "");
        BasepairState bpstate = str_to_basepairstate(bp_spl[1]);
        //Hack to stop memory out of bounds
        //TODO look into why this is happening!
        if(bp_spl.size() == 2) { bp_spl.push_back("c..."); }
        BasepairOP bp ( new Basepair(res1, res2, bpstate.r(), bp_spl[2] ));

        basepairs_.push_back(bp);
    }
    
    Strings end_indexes = split_str_by_delimiter(spl[7], " ");
    for (auto const & index : end_indexes) {
        ends_.push_back( basepairs_ [ std::stoi(index) ]);
    }
    end_ids_ = split_str_by_delimiter(spl[8], " ");
    secondary_structure_ = std::make_shared<sstruct::SecondaryStructure>(sstruct::str_to_secondary_structure(spl[9]));
}


Motif
Motif::copy() {
    Motif cmotif;
    cmotif.name_ = name_;
    cmotif.path_ = path_;
    cmotif.score_ = score_;
    cmotif.mtype_ = mtype_;
    cmotif.structure_ = StructureOP (new Structure(structure_->copy()));
    cmotif.beads_ = Beads(beads_.size());
    cmotif.basepairs_ = BasepairOPs();
    cmotif.ends_ = BasepairOPs();
    cmotif.secondary_structure_ = std::make_shared<sstruct::SecondaryStructure>(secondary_structure_->copy());
    cmotif.end_ids_ = end_ids_;
    cmotif.block_end_add_ = block_end_add_;
    int i = 0;
    for (auto const & b : beads_) {
        cmotif.beads_[i] = b.copy();
        i++;
    }
    i = 0;
    for (auto const & bp : basepairs_) {
        ResidueOP res1 = cmotif.get_residue(bp->res1()->uuid());
        ResidueOP res2 = cmotif.get_residue(bp->res2()->uuid());
        BasepairOP new_bp ( new Basepair ( res1, res2, bp->r(), bp->bp_type() )) ;
        new_bp->flipped(bp->flipped());
        new_bp->uuid(bp->uuid());
        cmotif.basepairs_.push_back(new_bp);
    }
    i = 0;
    for (auto & end: ends_) {
        BasepairOPs bps = cmotif.get_basepair(end->uuid());
        cmotif.ends_.push_back(bps[0]);
    }
    
    return cmotif;
}


String const
Motif::to_str() {
    std::stringstream ss;
    ss << path_ << "&" << name_ << "&" << score_ << "&" << block_end_add_ << "&" << mtype_;
    ss << "&" << structure_->to_str() << "&";
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
        ss << pos << " ";
    }
    ss << "&";
    for (auto const & end_id : end_ids_) {
        ss << end_id << " ";
    }
    ss << "&";
    ss << secondary_structure_->to_str();
    ss << "&";
    return ss.str();
}

String const
Motif::to_pdb_str() {
    return structure_->to_pdb_str();
}

void
Motif::to_pdb(String const fname) {
    return structure_->to_pdb(fname);
}

BasepairOPs
Motif::get_basepair(Uuid const & bp_uuid) {
    BasepairOPs bps;
    for( auto const & bp : basepairs_) {
        if(bp->uuid() == bp_uuid) { bps.push_back(bp); }
        if(bp->res1()->uuid() == bp_uuid || bp->res2()->uuid() == bp_uuid ) { bps.push_back(bp); }

    }
    return bps;
}

BasepairOPs
Motif::get_basepair(
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
Motif::get_basepair(
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
Motif::get_beads(
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
Motif::get_beads(
    BasepairOP const & excluded_end) {
    
    ResidueOPs excluded_res;
    excluded_res.push_back(excluded_end->res1());
    excluded_res.push_back(excluded_end->res2());
    beads_ = structure_->get_beads(excluded_res);
    return beads_;
}

BasepairOPs 
Motif::get_basepair(
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
Motif::end_index(BasepairOP const & end) {
    int pos = (int)(std::find(ends_.begin(), ends_.end(), end) - ends_.begin());
    return pos;
}


void
align_motif(
    BasepairStateOP const & ref_bp_state,
    BasepairOP const & motif_end,
    MotifOP & motif) {
    
    Matrix ref_T;
    transpose(ref_bp_state->r(), ref_T);
    Matrix r;
    dot(ref_T, motif_end->r(), r);
    Point trans = -motif_end->d();
    Transform t(r, trans);
    motif->transform(t);
    Point bp_pos_diff = ref_bp_state->d() - motif_end->d();
    motif->move(bp_pos_diff);
    bp_pos_diff = ref_bp_state->d() - motif_end->d();
    
    //align sugars for better overlap
    float dist1 = motif_end->res1()->get_atom("C1'")->coords().distance(ref_bp_state->sugars()[0]);
    float dist2 = motif_end->res2()->get_atom("C1'")->coords().distance(ref_bp_state->sugars()[1]);
    
    if (dist1 > 5 && dist2 > 5) { return; }
    
    Point sugar_diff_1, sugar_diff_2;
    if( dist1 < dist2 ) {
        sugar_diff_1 = ref_bp_state->sugars()[0] - motif_end->res1()->get_atom("C1'")->coords();
        sugar_diff_2 = ref_bp_state->sugars()[1]- motif_end->res2()->get_atom("C1'")->coords();
    }
    else {
        sugar_diff_1 = ref_bp_state->sugars()[0] - motif_end->res2()->get_atom("C1'")->coords();
        sugar_diff_2 = ref_bp_state->sugars()[1] - motif_end->res1()->get_atom("C1'")->coords();
    }
    
    motif->move( (sugar_diff_1 + sugar_diff_2) / 2);
}


MotifOP
get_aligned_motif(
    BasepairOP const & ref_bp,
    BasepairOP const & motif_end,
    MotifOP const & motif) {
    
    int motif_end_index = motif->end_index(motif_end);
    auto m_copy = std::make_shared<Motif>(motif->copy());
    auto new_motif_end = m_copy->ends()[motif_end_index];
    align_motif(ref_bp->state(), motif_end, m_copy);
    return m_copy;
    
}

MotifOP
file_to_motif(
    String const & path ) {
    Strings lines = get_lines_from_file(path);
    return std::make_shared<Motif>(lines[0],
                                   ResidueTypeSetManager::getInstance().residue_type_set());
}

Motif
ref_motif() {
    String path = resources_path() + "start.motif";
    String line;
    std::ifstream in;
    in.open(path);
    getline(in, line);
    in.close();
    if(line.size() == 0) {
        std::fstream infile;
        infile.open(path);
        infile >> line;
        std::cout << line << std::endl;
        exit(0);
    }
    Motif m ( line, ResidueTypeSetManager::getInstance().residue_type_set());
    return m;
}


















