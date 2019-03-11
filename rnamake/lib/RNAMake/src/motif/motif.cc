//
//  motif.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include <assert.h>

//RNAMake Headers
#include "base/file_io.h"
#include "base/settings.h"
#include "util/x3dna.h"
#include "structure/residue_type_set_manager.h"
#include "structure/chain.h"
#include "motif/motif.h"

Motif::Motif(
    String const & s,
    ResidueTypeSet const & rts):
RNAStructure(),
id_(Uuid()) {
 
    if(s.length() < 10) {
        throw "tried to construct Motif object from string, with a string too short";
    }
    
    auto spl       = base::split_str_by_delimiter(s, "&");
    path_          = spl[0];
    name_          = spl[1];
    score_         = std::stof(spl[2]);
    block_end_add_ = std::stoi(spl[3]);
    mtype_         = static_cast<MotifType>(std::stoi(spl[4]));
    structure_     = std::make_shared<Structure>(spl[5], rts);
    auto basepair_str = base::split_str_by_delimiter(spl[6], "@");
    for (auto const & bp_str : basepair_str) {
        auto bp_spl = base::split_str_by_delimiter(bp_str, ",");
        auto res_spl = base::split_str_by_delimiter(bp_spl[0], "-");
        auto res1_id = res_spl[0].substr(0,1);
        auto res2_id = res_spl[1].substr(0,1);
        auto res1_num = std::stoi(res_spl[0].substr(1));
        auto res2_num = std::stoi(res_spl[1].substr(1));
        auto res1 = structure_->get_residue(res1_num, res1_id, "");
        auto res2 = structure_->get_residue(res2_num, res2_id, "");
        auto bpstate = str_to_basepairstate(bp_spl[1]);
        //Hack to stop memory out of bounds
        //TODO look into why this is happening!
        if(bp_spl.size() == 2) { bp_spl.push_back("c..."); }
        BasepairOP bp ( new Basepair(res1, res2, bpstate.r(), bp_spl[2] ));

        basepairs_.push_back(bp);
    }
    
    Strings end_indexes = base::split_str_by_delimiter(spl[7], " ");
    for (auto const & index : end_indexes) {
        ends_.push_back( basepairs_ [ std::stoi(index) ]);
    }
    end_ids_ = base::split_str_by_delimiter(spl[8], " ");
    secondary_structure_ = std::make_shared<sstruct::Motif>(spl[9]);
    
    auto ss_res = secondary_structure_->residues();
    int i = 0;
    for(auto const & r : residues()) {
        ss_res[i]->uuid(r->uuid());
        i++;
    }
    
    if(spl.size() > 10) {
        auto bead_spl = base::split_str_by_delimiter(spl[10], ";");
        for(auto const & b_spl : bead_spl) {
            if(b_spl.length() < 5) { continue; }
            protein_beads_.push_back(Bead(b_spl));
        }
    }
}

Motif::Motif(
    Motif const & m) {
    
    name_       = m.name_;
    path_       = m.path_;
    score_      = m.score_;
    mtype_      = m.mtype_;
    structure_  = std::make_shared<Structure>(*m.structure_);
    beads_      = Beads(m.beads_.size());
    basepairs_  = BasepairOPs(m.basepairs_.size());
    ends_       = BasepairOPs(m.ends().size());
    end_ids_    = m.end_ids_;
    id_         = m.id_;
    block_end_add_ = m.block_end_add_;
    secondary_structure_ = std::make_shared<sstruct::Motif>(*m.secondary_structure_);
    protein_beads_ = Beads(m.protein_beads_.size());

    int i = 0;
    for (auto const & b : m.beads_) {
        beads_[i] = Bead(b);
        i++;
    }
    i = 0;
    for (auto const & b : m.protein_beads_) {
        protein_beads_[i] = Bead(b);
        i++;
    }
    i = 0;
    for (auto const & bp : m.basepairs_) {
        ResidueOP res1 = get_residue(bp->res1()->uuid());
        ResidueOP res2 = get_residue(bp->res2()->uuid());
        auto new_bp = std::make_shared<Basepair>(res1, res2, bp->r(), bp->bp_type());
        new_bp->flipped(bp->flipped());
        new_bp->uuid(bp->uuid());
        basepairs_[i] = new_bp;
        i++;
    }
    i = 0;
    for (auto & end : m.ends_) {
        auto bps = get_basepair(end->uuid());
        ends_[i] = bps[0];
        i++;
    }
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

void
Motif::new_res_uuids() {
    id_ = Uuid();
    for(auto & r : residues()) {
        auto ss_r = secondary_structure_->get_residue(r->uuid());
        r->new_uuid();
        ss_r->uuid(r->uuid());
    }
    for(auto & bp : basepairs()) { bp->uuid(Uuid()); }
}

void
Motif::copy_uuids_from_motif(
    Motif const & m) {
    
    id_ = m.id_;
    for(auto const & r : m.residues()) {
        auto this_r = get_residue(r->num(), r->chain_id(), r->i_code());
        this_r->uuid(r->uuid());
    }
    
    for(auto const & bp : m.basepairs()) {
        auto this_bp = get_basepair(bp->res1()->uuid(), bp->res2()->uuid())[0];
        this_bp->uuid(bp->uuid());
    }
}

MotifStateOP
Motif::get_state() {
    auto beads = get_beads(ends_[0]);
    Points bead_centers;
    for(auto const & b : beads) {
        if(b.btype() == BeadType::PHOS) { continue; }
        bead_centers.push_back(b.center());
    }
    BasepairStateOPs end_states;
    Strings end_names;
    for(auto const & end : ends_) {
        end_states.push_back(end->state());
        end_names.push_back(end->name());
    }
    MotifState ms(name_, end_names, end_ids_, end_states, bead_centers, score_,
                  (int)residues().size(), block_end_add_, id_);
    return std::make_shared<MotifState>(ms);
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
    r.unitarize();
    Point trans = -motif_end->d();
    Transform t(r, trans);
    motif->transform(t);
    Point bp_pos_diff = ref_bp_state->d() - motif_end->d();
    motif->move(bp_pos_diff);
    
    //align sugars for better overlap
    auto dist1 = motif_end->res1()->get_atom("C1'")->coords().distance(ref_bp_state->sugars()[0]);
    auto dist2 = motif_end->res2()->get_atom("C1'")->coords().distance(ref_bp_state->sugars()[0]);
    
    motif->get_beads(motif_end);
    if (dist1 > 5 && dist2 > 5) { return; }
    
    Point sugar_diff_1, sugar_diff_2;
    if( dist1 < dist2 ) {
        sugar_diff_1 = ref_bp_state->sugars()[0] - motif_end->res1()->get_atom("C1'")->coords();
        sugar_diff_2 = ref_bp_state->sugars()[1] - motif_end->res2()->get_atom("C1'")->coords();
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
    
    int motif_end_index = motif->get_end_index(motif_end);
    auto m_copy = std::make_shared<Motif>(*motif);
    auto new_motif_end = m_copy->ends()[motif_end_index];
    //align_motif(ref_bp->state(), motif_end, m_copy);
    align_motif(ref_bp->state(), new_motif_end, m_copy);
    return m_copy;
    
}

MotifOP
file_to_motif(
    String const & path ) {
    Strings lines = base::get_lines_from_file(path);
    return std::make_shared<Motif>(lines[0],
                                   ResidueTypeSetManager::getInstance().residue_type_set());
}

Motif
ref_motif() {
    String path = base::resources_path() + "start.motif";
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



















