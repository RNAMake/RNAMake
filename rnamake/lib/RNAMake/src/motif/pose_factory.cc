//
//  pose_factory.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 8/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>

#include "structure/chain.h"
#include "motif/pose_factory.h"
#include "motif/motif.h"

PoseOP
PoseFactory::pose_from_motif_tree(
    StructureOP const & structure,
    BasepairOPs const & basepairs,
    MotifOPs const & motifs,
    std::map<Uuid, int, UuidCompare> const & designable) {
    
    auto ends = mf_._setup_basepair_ends(structure, basepairs);
    auto m = std::make_shared<Motif>(structure, basepairs, ends);
    mf_._setup_secondary_structure(m);
    
    auto p = std::make_shared<Pose>(m);
    _add_motifs_to_pose(p, motifs);

    return p;
    
}


void
PoseFactory::_add_motifs_to_pose(
    PoseOP & p,
    MotifOPs const & motifs) {
    
    //int j = 0;
    float dist, best_dist;
    Point r1_cent, r2_cent;
    ResidueOP best_match;
    for(auto const & m : motifs) {
        BasepairOPs bps;
        std::map<ResidueOP, ResidueOP> residue_map;
        
        for(auto const & bp : m->basepairs()) {
            auto new_bp = p->get_basepair(bp->uuid());
            if(new_bp.size() != 0) {
                bps.push_back(new_bp[0]);
                continue;
            }
            
            int best = 1000;
            BasepairOP best_bp = nullptr;
            for(auto const & m2 : motifs) {
                if(m == m2) { continue; }
                for(auto const & end : m2->ends()) {
                    auto alt_bp = p->get_basepair(end->uuid());
                    if(alt_bp.size() == 0) { continue; }
                    dist = bp->d().distance(alt_bp[0]->d());
                    if(dist < best) {
                        best_bp = alt_bp[0];
                        best = dist;
                    }
                }
            }
            
            if(best_bp == nullptr) { continue; }
            new_bp = BasepairOPs { best_bp };
            for(auto const & r1 : best_bp->residues()) {
                r1_cent = center(r1->atoms());
                best_dist = 1000;
                best_match = nullptr;
                for (auto const & r2 : bp->residues()) {
                    r2_cent = center(r2->atoms());
                    dist = r1_cent.distance(r2_cent);
                    if(dist < best_dist) {
                        best_dist = dist;
                        best_match = r2;
                    }
                }
                residue_map[best_match] = r1;
            }
            bps.push_back(new_bp[0]);
        }
        
        ChainOPs chains;
        for(auto const & c : m->chains()) {
            ResidueOPs res;
            for(auto const & r : c->residues()) {
                auto new_r = p->get_residue(r->uuid());
                if     (new_r != nullptr) {
                    res.push_back(new_r);
                }
                else if(residue_map.find(r) != residue_map.end()) {
                    res.push_back(residue_map[r]);
                }
                else {
                    throw std::runtime_error("cannot find residue in PoseFactory::_add_motifs_to_pose");
                }
            }
            chains.push_back(std::make_shared<Chain>(res));
        }
    }

    
    
}




































