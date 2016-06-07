//
//  sample_helix.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 6/7/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "sample_helix.hpp"

#include "base/cl_option.h"
#include "secondary_structure/ss_tree.h"
#include "motif/motif_tree.h"
#include "resources/resource_manager.h"
#include "motif/motif_factory.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"


int main(int argc, const char * argv[]) {
    auto seq = String("GGCCCUCAAGGG&CCCUUGAGGGCC");
    auto ss  = String("((((((((((((&))))))))))))");
    
    sstruct::SS_Tree ss_tree(seq, ss);
    sstruct::SS_TreeNodeOP last_node = nullptr;

    auto mt = std::make_shared<MotifTree>();
    
    int i = -1;
    for(auto const & n : ss_tree) {
        i++;
        if(i == 0) { continue; }
        if(n->data()->sequence() == "&&&") { continue; }
        if(last_node == nullptr ) { last_node = n; continue; }
        
        auto seq1 = last_node->data()->sequence();
        auto seq2 = n->data()->sequence();
        auto motif_name = String("");
        motif_name.push_back(seq1[0]); motif_name.push_back(seq1[2]);
        motif_name.push_back('=');
        motif_name.push_back(seq2[0]); motif_name.push_back(seq2[2]);

        last_node = n;
        mt->add_motif(ResourceManager::getInstance().get_motif(motif_name));
        
    }
    
    auto mf = MotifFactory();
    auto m = mf.motif_from_file("/Users/josephyesselman/Downloads/sele.pdb");
    m->block_end_add(-1);
    m->ends()[0]->flip();
    m->mtype(HAIRPIN);
    ResourceManager::getInstance().register_motif(m);
    mt->add_motif(m);

    auto mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);
    
    auto sampler = ThermoFlucSampler();
    sampler.setup(mset);
    sampler.to_pdb("start.pdb");
    
    for(int i = 0; i < 10000000; i++) {
        try {
            sampler.next();
        } catch(...) { }
            
        if(i % 10000 == 0 ) {
            try {
                sampler.to_pdb("sampled."+std::to_string(i)+".pdb");
            } catch(...) { }
        }
        
    }

    
    return 0;
}