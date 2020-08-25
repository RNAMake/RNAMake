//
//  motif_factory.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>

#include "motif/motif_factory.h"
#include "structure/chain.h"
#include "secondary_structure/util.h"
#include "base/file_io.h"
#include "util/x3dna.h"


namespace motif {

MotifOP
MotifFactory::motif_from_file(
        String const & path,
        bool rebuild_x3dna,
        bool include_protein,
        int force_num_chains) {


    if (!base::file_exists(path)) {
        throw MotifFactoryException("cannot generate motif from file " + path + " does not exist");
    }
    parser_ = MotiftoSecondaryStructure();
    auto fname = base::filename(path);
    auto pdb_path = path;
    structure::StructureOP structure;
    if (base::is_dir(path)) {
        rebuild_x3dna = false;

        pdb_path = path + "/" + fname + ".pdb";
        if (!base::file_exists(pdb_path)) {
            throw MotifFactoryException(
                    "cannot generate motif from directory " + path + " it exists but there is no pdb "
                            " in it expected: dir_name/dir_name.pdb");
        }

        structure = std::make_shared<structure::Structure>(pdb_path);
    } else {
        structure = std::make_shared<structure::Structure>(path);
        fname = fname.substr(0, fname.length() - 4);
    }

    
    if (force_num_chains != -1 && structure->chains().size() != force_num_chains) {
        if (force_num_chains > structure->chains().size()) {
            throw MotifFactoryException(
                    "force_num_chains must be smaller than current number in motif");
        }
        structure = _get_reduced_chain_num_structure(*structure, force_num_chains);
    }

    auto basepairs = _setup_basepairs(pdb_path, structure, rebuild_x3dna);
    auto ends = _setup_basepair_ends(structure, basepairs);
    auto m = std::make_shared<Motif>(structure, basepairs, ends);
    m->name(fname);
    m->path(path);

    //clean up x3dna generated files if they were created
    try {
        std::remove("ref_frames.dat");
        std::remove(String(fname + "_dssr.out").c_str());
    } catch (...) {}
    _setup_secondary_structure(m);

    if (include_protein) {
        auto res = pdb_parser_.parse(pdb_path, true, false);
        auto beads = structure::Beads();
        for (auto const & r : res) {
            try {
                auto bead = structure::Bead(r->get_atom("CA")->coords(), structure::BeadType::BASE);
                beads.push_back(bead);
            } catch (...) { continue; }
        }
        m->protein_beads(beads);

    }

    if (m->ends().size() > 0) {
        return get_oriented_motif(m, 0);
    } else {
        return m;
    }
}


MotifOP
MotifFactory::get_oriented_motif(
        MotifOP const & m,
        int end_pos) {

    auto m_added = can_align_motif_to_end(m, end_pos);
    m_added = align_motif_to_common_frame(m_added, end_pos);
    _standardize_motif(m_added);
    auto m_added_2_1 = get_aligned_motif(m->ends()[end_pos], m_added->ends()[0], m_added);

    auto count_1 = _bead_overlap(m, m_added_2_1);
    if (count_1 == 0) { return m_added_2_1; }

    m->ends()[end_pos]->flip();
    auto m_added_2_2 = get_aligned_motif(m->ends()[end_pos], m_added->ends()[0], m_added);
    auto count_2 = _bead_overlap(m, m_added_2_2);

    if (count_2 == 0) { return m_added_2_2; }

    if (count_1 > count_2) { return m_added_2_2; }
    else { return m_added_2_1; }

}

MotifOP
MotifFactory::motif_from_res(
        structure::ResidueOPs & res,
        structure::BasepairOPs const & bps) {

    auto chains = structure::ChainOPs();
    connect_residues_into_chains(res, chains);
    auto structure = std::make_shared<structure::Structure>(chains);
    auto ends = _setup_basepair_ends(structure, bps);

    if (ends.size() == 0 && bps.size() >= 2) {
        std::cout << ends.size() << " " << bps.size() << " " << chains.size() << std::endl;
        structure->to_pdb("test.pdb");
        throw MotifFactoryException(
                "unexpected number of ends when generating a motif from residues");
    }

    auto m = std::make_shared<Motif>(structure, bps, ends);
    _setup_secondary_structure(m);
    return m;
}

MotifOP
MotifFactory::motif_from_bps(
        structure::BasepairOPs const & bps) {

    auto res = structure::ResidueOPs();
    for (auto const & bp : bps) {
        res.push_back(bp->res1());

        res.push_back(bp->res2());
    }
    auto chains = structure::ChainOPs();
    connect_residues_into_chains(res, chains);
    auto structure = std::make_shared<structure::Structure>(chains);
    auto ends = _setup_basepair_ends(structure, bps);

    if (ends.size() != 2 && bps.size() >= 2) {
        throw MotifFactoryException(
                "unexpected number of ends when generating a motif from basepairs");
    }

    auto m = std::make_shared<Motif>(structure, bps, ends);
    _setup_secondary_structure(m);
    return m;

}


MotifOP
MotifFactory::can_align_motif_to_end(
        MotifOP const & m,
        int ei) {

    auto m_end = m->ends()[ei];
    auto m_updated = std::make_shared<Motif>(*m);
    auto m_added = MotifOP(nullptr);
    auto m_added_1 = get_aligned_motif(base_motif_->ends()[1], m_end, m);
    m_end->flip();
    auto m_added_2 = get_aligned_motif(base_motif_->ends()[1], m_end, m);

    auto clash_count_1 = _steric_clash_count(base_motif_, m_added_1);
    auto clash_count_2 = _steric_clash_count(base_motif_, m_added_2);

    if (clash_count_1 > clash_count_2) {
        m_updated->ends()[ei]->flip();
        m_added = m_added_2;
    } else {
        m_added = m_added_1;
    }

    auto i = -1;
    for (auto & end : m_added->ends()) {
        i++;
        if (end->name() == m_end->name()) { continue; }
        auto added_helix_1 = get_aligned_motif(end, added_helix_->ends()[0], added_helix_);
        end->flip();
        auto added_helix_2 = get_aligned_motif(end, added_helix_->ends()[0], added_helix_);
        end->flip();

        clash_count_1 = _steric_clash_count(m_added, added_helix_1) + _steric_clash_count(base_motif_, added_helix_1);
        clash_count_2 = _steric_clash_count(m_added, added_helix_2) + _steric_clash_count(base_motif_, added_helix_2);

        if (clash_count_1 > clash_count_2) {
            m_updated->ends()[i]->flip();
        }
    }

    return m_updated;
}


void
MotifFactory::standardize_rna_structure_ends(
        MotifOP & m) {

    m->get_beads(m->ends());
    for (auto const & end : m->ends()) {
        auto m_added = get_aligned_motif(end, added_helix_->ends()[0], added_helix_);
        if (_steric_clash(m, m_added) == 0) { continue; }
        else { end->flip(); }
    }

}

MotifOP
MotifFactory::align_motif_to_common_frame(
        MotifOP const & m,
        int ei) {

    auto m_added = get_aligned_motif(ref_motif_->ends()[0], m->ends()[ei], m);
    return m_added;
}

void
MotifFactory::_standardize_motif(
        MotifOP & m) {

    _align_chains(m);
    _align_ends(m);
    _setup_secondary_structure(m);
}


structure::BasepairOPs
MotifFactory::_setup_basepairs(
        String const & path,
        structure::StructureOP const & s,
        bool rebuild_x3dna) {

    auto basepairs = structure::BasepairOPs();
    auto x3dna_parser = util::X3dna();

#ifdef JSON_BASEPAIRS 
    auto x_basepairs = x3dna_parser.get_basepairs_json(path);
    structure::ResidueOP res1, res2;
    //BasepairOP bp;
    String i_code_1, i_code_2;

    for (auto const & xbp : x_basepairs) {
        // super hacky way of converting back to the old string system.
        if(xbp.bp_type != util::X3dnaBPType::cWUW) {
            continue;
        } 
        
        if(xbp.res1.i_code == ' ') {
            i_code_1 = "";
        }
        else {
            i_code_1 = String(1, xbp.res1.i_code);
        }

        if(xbp.res2.i_code == ' ') {
            i_code_2 = "";
        }
        else {
            i_code_2 = String(1, xbp.res2.i_code);
        }
        

        res1 = s->get_residue(xbp.res1.num, String(1, xbp.res1.chain_id), i_code_1);
        res2 = s->get_residue(xbp.res2.num, String(1, xbp.res2.chain_id), i_code_2);
        if (res1 == nullptr || res2 == nullptr) {
            //std::cout << xbp.res1.num << " " << xbp.res2.num << std::endl;
            throw MotifFactoryException("cannot find residues in basepair during setup");
        }

        // Emplacement constructs the object in-place
        basepairs.emplace_back(new structure::Basepair(res1, res2, xbp.r, util::get_str_from_x3dna_type(xbp.bp_type)));
    }

#else 
    auto x_basepairs = x3dna_parser.get_basepairs(path);
    structure::ResidueOP res1, res2;
    //BasepairOP bp;
    String i_code_1, i_code_2;

    for (auto const & xbp : x_basepairs) {
        // super hacky way of converting back to the old string system.
        if(xbp.res1.i_code == ' ') {
            i_code_1 = "";
        }
        else {
            i_code_1 = String(1, xbp.res1.i_code);
        }

        if(xbp.res2.i_code == ' ') {
            i_code_2 = "";
        }
        else {
            i_code_2 = String(1, xbp.res2.i_code);
        }

        res1 = s->get_residue(xbp.res1.num, String(1, xbp.res1.chain_id), i_code_1);
        res2 = s->get_residue(xbp.res2.num, String(1, xbp.res2.chain_id), i_code_2);
        if (res1 == nullptr || res2 == nullptr) {
            throw MotifFactoryException("cannot find residues in basepair during setup");
        }

        // Emplacement constructs the object in-place
        basepairs.emplace_back(new structure::Basepair(res1, res2, xbp.r, util::get_str_from_x3dna_type(xbp.bp_type)));
    }
#endif

    return basepairs;

}

structure::BasepairOPs
MotifFactory::_setup_basepair_ends(
        structure::StructureOP const & structure,
        structure::BasepairOPs const & basepairs) {

    structure::ResidueOPs chain_ends;
    for (auto const & c : structure->chains()) {
        chain_ends.push_back(c->first());
        if (c->residues().size() > 1) { chain_ends.push_back(c->last()); }
    }

    structure::BasepairOPs ends;
    for (auto const & bp : basepairs) {
        for (auto const & ce1 : chain_ends) {
            for (auto const & ce2 : chain_ends) {
                if (bp->bp_type() == "cW-W" &&
                    bp->res1()->uuid() == ce1->uuid() &&
                    bp->res2()->uuid() == ce2->uuid()) {
                    ends.push_back(bp);
                }
            }
        }
    }

    return ends;
}

void
MotifFactory::_setup_secondary_structure(
        MotifOP & m) {
    
    auto ss = parser_.to_secondary_structure(m);
    Strings end_ids(m->ends().size());
    int i = 0;
    for (auto const & end : m->ends()) {
        auto res1 = ss->get_residue(end->res1()->uuid());
        auto res2 = ss->get_residue(end->res2()->uuid());
        auto ss_end = ss->get_basepair(res1, res2)[0];
        end_ids[i] = secondary_structure::assign_end_id(ss, ss_end);
        i++;
    }

    ss->end_ids(end_ids);
    m->end_ids(end_ids);
    m->secondary_structure(std::static_pointer_cast<secondary_structure::Motif>(ss));

}

void
MotifFactory::_align_chains(
        MotifOP & m) {

    auto chains = m->chains();
    structure::ChainOP closest;
    float best = 1000;
    auto c2 = structure::center(ref_motif_->ends()[0]->atoms());
    math::Point c1;
    for (auto const & c : chains) {
        c1 = structure::center(c->first()->atoms());
        float dist = c1.distance(c2);
        if (dist < best) {
            best = dist;
            closest = c;
        }
    }

    auto updated_chains = structure::ChainOPs{closest};
    for (auto const & c : chains) {
        if (c != closest) { updated_chains.push_back(c); }
    }

    m->structure(std::make_shared<structure::Structure>(updated_chains));
}

void
MotifFactory::_align_ends(
        MotifOP & m) {

    structure::BasepairOP closest;
    float best = 1000;
    auto c2 = structure::center(ref_motif_->ends()[0]->atoms());
    math::Point c1;

    for (auto const & end : m->ends()) {
        c1 = structure::center(end->atoms());
        float dist = c1.distance(c2);
        if (dist < best) {
            best = dist;
            closest = end;
        }
    }

    auto updated_ends = structure::BasepairOPs{closest};
    for (auto const & end : m->ends()) {
        if (end != closest) { updated_ends.push_back(end); }
    }

    for (auto & end : m->ends()) {
        int flip_res = 0;
        for (auto const & c : m->chains()) {
            if (c->first() == end->res2()) {
                flip_res = 1;
                break;
            }
        }

        if (!flip_res) { continue; }

        auto temp = end->res1();
        end->res1(end->res2());
        end->res2(temp);
    }

    m->ends(updated_ends);

}

int
MotifFactory::_steric_clash(
        MotifOP const & m1,
        MotifOP const & m2) {

    for (auto const & c1 : m1->beads()) {
        for (auto const & c2 : m2->beads()) {
            if (c1.btype() == structure::BeadType::PHOS || c2.btype() == structure::BeadType::PHOS) { continue; }
            float dist = c1.center().distance(c2.center());
            if (dist < clash_radius_) { return 1; }
        }
    }

    return 0;

}

int
MotifFactory::_steric_clash_count(
        MotifOP const & m1,
        MotifOP const & m2) {

    auto count = 0;
    for (auto const & c1 : m1->beads()) {
        for (auto const & c2 : m2->beads()) {
            if (c1.btype() == structure::BeadType::PHOS || c2.btype() == structure::BeadType::PHOS) { continue; }
            float dist = c1.center().distance(c2.center());
            if (dist < clash_radius_) { count += 1; }
        }
    }

    return count;

}

int
MotifFactory::_bead_overlap(
        MotifOP const & m1,
        MotifOP const & m2) {

    auto count = 0;
    auto dist = 0.0;
    for (int i = 0; i < m1->beads().size(); i++) {
        dist = m1->beads()[i].distance(m2->beads()[i]);
        if (dist < 0.2) { count += 1; }
    }

    return count;

}

structure::StructureOP
MotifFactory::_get_reduced_chain_num_structure(
        structure::Structure const & start,
        int chain_num) {

    auto chains = start.chains();

    for (int round = 0; round < (chains.size() - chain_num); round++) {
        auto best_i = -1;
        auto best_j = -1;
        auto best_score = 1000;
        auto dist = 100;

        for (int i = 0; i < chains.size(); i++) {
            for (int j = 0; j < chains.size(); j++) {
                if (i == j) { continue; }
                auto o3_atom = chains[i]->last()->get_atom("O3'");
                auto p_atom = chains[j]->first()->get_atom("P");
                dist = o3_atom->coords().distance(p_atom->coords());
                if (dist < best_score) {
                    best_i = i;
                    best_j = j;
                    best_score = dist;
                }

            }
        }

        auto new_chains = structure::ChainOPs();
        auto residues = structure::ResidueOPs();
        for (auto const & r : chains[best_i]->residues()) {
            residues.push_back(r);
        }
        for (auto const & r : chains[best_j]->residues()) {
            residues.push_back(r);
        }
        auto new_chain = std::make_shared<structure::Chain>(residues);
        new_chains.push_back(new_chain);
        for (int i = 0; i < chains.size(); i++) {
            if (i == best_i || i == best_j) { continue; }
            new_chains.push_back(chains[i]);
        }
        chains = new_chains;
    }

    return std::make_shared<structure::Structure>(chains);
}

}

















