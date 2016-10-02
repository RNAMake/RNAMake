//
//  sequence_optimizer_3d.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/30/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "util/monte_carlo.h"
#include "motif_data_structures/motif_state_tree.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"


SequenceOptimizer3D::SequenceOptimizer3D():
sequence_designer_(eternabot::SequenceDesigner()),
rng_(RandomNumberGenerator())
{}


SequenceOptimizer3D::OptimizedSequenceOPs
SequenceOptimizer3D::get_optimized_sequences(
    MotifTreeOP const & mt,
    BasepairOP const & target_bp,
    int ni,
    int ei) {
    
    auto sols = OptimizedSequenceOPs();
    auto ss = mt->designable_secondary_structure();
    
    auto designable_bps = DesignableBPOPs();
    for(auto const & bp : ss->basepairs()) {
        auto bp_name = bp->res1()->name()+bp->res2()->name();
        if(bp_name == "NN") {
            designable_bps.push_back(
                std::make_shared<DesignableBP>(DesignableBP{bp, "", nullptr, nullptr}));
        }
    }
    
    auto possible_bps = std::vector<Strings>({{"A", "U"}, {"U", "A"}, {"G", "C"}, {"C", "G"}});
    
    for(auto const & d_bp : designable_bps) {
        for(auto const & m : ss->motifs()) {
            if(m->name() != "HELIX.IDEAL") { continue; }
            if(m->ends()[0] == d_bp->bp) { d_bp->m_id_bot = std::make_shared<Uuid>(m->id()); }
            if(m->ends()[1] == d_bp->bp) { d_bp->m_id_top = std::make_shared<Uuid>(m->id()); }
        }
        
        auto state = possible_bps[rng_.randrange(possible_bps.size())];
        d_bp->bp->res1()->name(state[0]);
        d_bp->bp->res2()->name(state[1]);
    }
    
    auto mst = std::make_shared<MotifStateTree>(mt);
    for(auto const & m : ss->motifs()) {
        if(m->name() != "HELIX.IDEAL") { continue; }
        auto n = mst->get_node(m->id());
    }
    
    
    
    
    return sols;
    
}