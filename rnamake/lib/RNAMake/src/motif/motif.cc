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
#include "util/settings.h"
#include "util/x3dna.h"
#include "structure/resource_manager.h"
#include "structure/chain.h"
#include "motif/motif.h"

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
 
    if(s.length() < 10) {
        throw "tried to construct Motif object from string, with a string too short";
    }
    
    Strings spl = split_str_by_delimiter(s, "&");
    mdir_ = spl[0];
    name_ = spl[1];
    score_ = std::stof(spl[2]);
    mtype_ = static_cast<MotifType>(std::stoi(spl[3]));
    structure_ = StructureOP( new Structure(str_to_structure(spl[4], rts)));
    Strings basepair_str = split_str_by_delimiter(spl[5], "@");
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
        BasepairOP bp ( new Basepair(res1, res2, bpstate.r(), bp_spl[2] ));
        bp->flipped(std::stoi(bp_spl[4]));

        basepairs_.push_back(bp);
    }
    
    Strings end_indexes = split_str_by_delimiter(spl[6], " ");
    for (auto const & index : end_indexes) {
        ends_.push_back( basepairs_ [ std::stoi(index) ]);
    }
    _cache_basepair_frames();
}

Motif::Motif(
    String const & path):
    beads_(Beads()),
    score_(0),
    basepairs_(BasepairOPs()),
    ends_(BasepairOPs()),
    mdir_(String()),
    name_(String()),
    cached_rotations_(Matrices()) {
        
    struct stat s;
    if( stat(path.c_str(),&s) == 0 ) {
        //it's a directory
        if( s.st_mode & S_IFDIR ) {
            String fname = filename(path);
            mdir_ = path;
            name_ = fname;
            structure_ = StructureOP(new Structure(mdir_ + "/" + fname + ".pdb"));
        }
        //it's a file
        else if( s.st_mode & S_IFREG )
        {
            mdir_ = base_dir(path);
            String fname = filename(path);
            name_ = fname.substr(0, fname.length()-4);
            structure_ = StructureOP(new Structure(path));
        }

    }

    _setup_basepairs();
    setup_basepair_ends();
    _cache_basepair_frames();

}



Motif
Motif::copy() {
    Motif cmotif;
    cmotif.name_ = name_;
    cmotif.mdir_ = mdir_;
    cmotif.score_ = score_;
    cmotif.mtype_ = mtype_;
    cmotif.structure_ = StructureOP (new Structure(structure_->copy()));
    cmotif.beads_ = Beads(beads_.size());
    cmotif.basepairs_ = BasepairOPs();
    cmotif.ends_ = BasepairOPs();
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
    
    cmotif._cache_basepair_frames();
    return cmotif;
}

void
Motif::setup_basepair_ends() {
    //TODO check to see if this works the same as python, its different
    
    ResidueOPs chain_ends;
    for(auto const & c : chains()) {
        chain_ends.push_back(c->first());
        if(c->residues().size() > 1) { chain_ends.push_back(c->last()); }
    }
    
    for( auto const & bp : basepairs_ ) {
        for (auto const & ce1 : chain_ends) {
            for(auto const & ce2 : chain_ends) {
                if(bp->bp_type().compare("cW-W") == 0 && bp->res1() == ce1 && bp->res2() == ce2) {
                    ends_.push_back(bp);
                }
                
            }
        }
    }
}

String const
Motif::to_str() {
    std::stringstream ss;
    ss << mdir_ << "&" << name_ << "&" << score_ << "&" << mtype_ << "&" << structure_->to_str() << "&";
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

Beads const &
Motif::get_beads(BasepairOPs const & excluded_ends) {
    ResidueOPs excluded_res;
    for (auto const & end : excluded_ends) {
        excluded_res.push_back(end->res1());
        excluded_res.push_back(end->res2());
    }
    beads_ = structure_->get_beads(excluded_res);
    return beads_;
}

Beads const &
Motif::get_beads(BasepairOP const & excluded_end) {
    ResidueOPs excluded_res;
    excluded_res.push_back(excluded_end->res1());
    excluded_res.push_back(excluded_end->res2());
    beads_ = structure_->get_beads(excluded_res);
    return beads_;
}

String
Motif::sequence() {
    String seq;
    ChainOPs const cs = chains();
    int i = -1;
    for (auto const & c : cs) {
        i++;
        for (auto const & r : c->residues()) {
            seq += r->short_name();
        }
        if(i+1 != cs.size()) {
            seq += "&";
        }
    }
    return seq;
}

String
Motif::secondary_structure() {
    String structure, ss;
    BasepairOPs bps;
    BasepairOP saved_bp(NULL);
    ResidueOP partner_res;
    std::map<String, int> seen_bp, seen_res;
    int is_bp = 0, passes = 0, bp_seen = 0, res_seen = 0, partner_res_seen = 0;
    int count = -1;
    for (auto const & c : chains()) {
        for (auto const & r : c->residues()) {
            count ++;
            ss = "";
            is_bp = 0;
            bps = get_basepair(r->uuid());
            for(auto const & bp : bps) {
                partner_res = bp->partner(r);
                is_bp = 1;
                passes = 0;
                saved_bp.reset();
                if(wc_bp(bp) && bp->bp_type().compare("cW-W") == 0) { passes = 1; }
                if(gu_bp(bp) && bp->bp_type().compare("cW-W") == 0) { passes = 1; }
                
                if(passes) {
                    saved_bp = bp;
                    bp_seen = 0; partner_res_seen = 0; res_seen = 0;
                    if(seen_bp.find(bp->uuid().s_uuid()) != seen_bp.end())  { bp_seen = 1;  }
                    if(seen_res.find(partner_res->uuid().s_uuid()) != seen_res.end()) { partner_res_seen = 1; }
                    if(seen_res.find(r->uuid().s_uuid()) != seen_res.end()) { res_seen = 1; }
                    if(bp_seen == 0 && res_seen == 0 && partner_res_seen == 0) {
                        seen_res[r->uuid().s_uuid()] = 1;
                        ss = "(";
                    }
                    else if(partner_res_seen == 1) {
                        if(seen_res[partner_res->uuid().s_uuid()] > 1) { ss = "."; }
                        else {
                            ss = ")";
                            seen_res[r->uuid().s_uuid()] = 1;
                            seen_res[partner_res->uuid().s_uuid()] += 1;
                            break;
                        }
                    }
                }
                else if(seen_res.find(r->uuid().s_uuid()) == seen_res.end()) {  ss = "."; }

            }
            if(!is_bp) { ss = "."; }
            if(saved_bp.get() != NULL) { seen_bp[saved_bp->uuid().s_uuid()] = 1; }
            structure += ss;
        }
        structure += "&";
    }
    
    return structure.substr(0, structure.length()-1);

}

void
Motif::_setup_basepairs() {
    
    X3dna x3dna_parser;
    String mdir = mdir_;
    if(mdir.length() < 3) { mdir = ""; }
    
    X3Basepairs x_basepairs = x3dna_parser.get_basepairs(mdir + "/" + name_);
    ResidueOP res1, res2;
    BasepairOP bp;
    for(auto const & xbp : x_basepairs) {
        res1 = structure_->get_residue(xbp.res1.num, xbp.res1.chain_id, xbp.res1.i_code);
        res2 = structure_->get_residue(xbp.res2.num, xbp.res2.chain_id, xbp.res2.i_code);
        if (res1 == nullptr || res2 == nullptr) {
            throw "cannot find residues in basepair during setup";
        }
        
        bp = BasepairOP(new Basepair(res1, res2, xbp.r, xbp.bp_type));
        _assign_bp_primes(bp);
        basepairs_.push_back(bp);
    }
}

void
Motif::_assign_bp_primes(BasepairOP & bp) {
    int res1_pos = 0, res2_pos = 0, res1_total = 0, res2_total = 0;
    int i = 0;
    for(auto const & c : chains()) {
        i = 0;
        for(auto const & r : c->residues()) {
            if(bp->res1() == r) {
                res1_pos = i;
                res1_total = (int)c->residues().size();
            }
            else if(bp->res2() == r) {
                res2_pos = i;
                res2_total = (int)c->residues().size();
            }
            i++;
        }
    }
    
    if(res1_pos > res1_total/2 || res2_pos < res2_total/2) {
        ResidueOP temp = bp->res1();
        bp->res1(bp->res2());
        bp->res2(temp);
    }
    
}

BasepairOP const &
Motif::get_basepair_by_name(
    String const & name) {
    Strings name_spl = split_str_by_delimiter(name, "-");
    String alt_name = name_spl[1] + "-" + name_spl[0];
    for(auto const & bp : basepairs_) {
        if(name.compare(bp->name()) == 0 || alt_name.compare(bp->name()) == 0) {
            return bp;
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
align_motif(BasepairOP const & ref_bp,
            BasepairOP const & motif_end,
            MotifOP const & motif) {
    
    Matrix ref_T;
    transpose(ref_bp->r(), ref_T);
    Matrix r;
    dot(ref_T, motif_end->r(), r);
    Point trans = -motif_end->d();
    Transform t(r, trans);
    motif->transform(t);
    Point bp_pos_diff = ref_bp->d() - motif_end->d();
    motif->move(bp_pos_diff);
    bp_pos_diff = ref_bp->d() - motif_end->d();
    
    //align sugars for better overlap
    float dist1 = motif_end->res1()->get_atom("C1'")->coords().distance(ref_bp->res1()->get_atom("C1'")->coords());
    float dist2 = motif_end->res2()->get_atom("C1'")->coords().distance(ref_bp->res1()->get_atom("C1'")->coords());
    
    if (dist1 > 5 && dist2 > 5) { return; }
    
    Point sugar_diff_1, sugar_diff_2;
    if( dist1 < dist2 ) {
        sugar_diff_1 = ref_bp->res1()->get_atom("C1'")->coords() - motif_end->res1()->get_atom("C1'")->coords();
        sugar_diff_2 = ref_bp->res2()->get_atom("C1'")->coords() - motif_end->res2()->get_atom("C1'")->coords();
    }
    else {
        sugar_diff_1 = ref_bp->res1()->get_atom("C1'")->coords() - motif_end->res2()->get_atom("C1'")->coords();
        sugar_diff_2 = ref_bp->res2()->get_atom("C1'")->coords() - motif_end->res1()->get_atom("C1'")->coords();
    }
    
    motif->move( (sugar_diff_1 + sugar_diff_2) / 2);
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
    Motif m ( line, ResourceManager::getInstance().residue_type_set());
    return m;
}

















