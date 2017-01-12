//
//  pose.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/pose.h"
#include "secondary_structure/util.h"

namespace sstruct {

void
Pose::_build_helices() {
    auto steps = MotifOPs();
    for(auto const & m : motifs_) {
        if(m->mtype() != MotifType::HELIX) { continue; }
        steps.push_back(m);
    }
    
    auto seen = std::map<MotifOP, int>();
    auto current = MotifOP(nullptr);
    auto helix_motifs = MotifOPs();
    int found = 0;
    while(1) {
        helix_motifs = MotifOPs();
        current = nullptr;
        for(auto const & m1 : steps) {
            found = 0;
            if(seen.find(m1) != seen.end()) { continue; }
            for(auto const & m2 : steps) {
                if(seen.find(m2) != seen.end()) { continue; }
                if(m1 == m2) { continue; }
                for(auto const & end : m2->ends()) {
                    if(m1->ends()[0] == end) {
                        found = 1;
                        break;
                    }
                }
            }
            if(!found) { current = m1; break; }
        }
        
        
        if(found) { break; }
        if(current == nullptr) { break; }
        
        
        found = 1;
        while(found) {
            seen[current] = 1;
            helix_motifs.push_back(current);
            found = 0;
            for(auto const & m : steps) {
                if(seen.find(m) != seen.end()) { continue; }
                if(m->ends()[0] == current->ends()[1]) {
                    current = m;
                    found = 1;
                    break;
                }
            }
            
        }
        
        ResidueOPs res1, res2;
        auto bps = BasepairOPs();
        auto ends = BasepairOPs();
        
        res1.push_back(helix_motifs[0]->chains()[0]->first());
        res2.push_back(helix_motifs[0]->chains()[1]->last());
        bps.push_back(helix_motifs[0]->ends()[0]);
        ends.push_back(helix_motifs[0]->ends()[0]);
        ends.push_back(helix_motifs.back()->ends()[1]);

 
        for(auto const & m : helix_motifs) {
            res1.push_back(m->chains()[0]->last());
            res2.push_back(m->chains()[1]->first());
            bps.push_back(m->ends()[1]);
        }
        
        std::reverse(res2.begin(), res2.end());
        
        auto chains = ChainOPs{std::make_shared<Chain>(res1),
                               std::make_shared<Chain>(res2) };
        auto struc = std::make_shared<Structure>(chains);
        auto h_m = std::make_shared<Motif>(struc, bps, ends);
        
        helices_.push_back(h_m);
        
    }
    
}

void
Pose::replace_sequence(
    String const & seq) {
    
    RNAStructure::replace_sequence(seq);

    for(auto & m : motifs_) {
        int i = 0;
        auto end_ids = Strings(m->ends().size());
        for(auto const & end : m->ends()) {
            end_ids[i] = assign_end_id(m, end);
            i++;
        }
        m->end_ids(end_ids);
    }
}

void
Pose::update_motif(Uuid const & uuid) {
    auto m = motif(uuid);
    auto end_ids = Strings(m->ends().size());
    int i = 0;
    for(auto const & end : m->ends()) {
        end_ids[i] = assign_end_id(m, end);
        i++;
    }
    m->end_ids(end_ids);
}

    
}

