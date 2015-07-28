//
//  motif_tree_topology_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_topology_unittest.h"
#include "util/random_number_generator.h"
#include "motif/motif_tree_topology.h"
#include "motif/motif_tree.h"
#include "resources/library_manager.h"

Strings
random_helix() {
    Strings pairs = split_str_by_delimiter("AU,UA,GC,CG",",");
    
    String s1,s2;
    String ss1, ss2;
    RandomNumberGenerator rng;
    //int size = rng.randrange(50);
    int size = 3;
    for(int i = 0; i < size; i++) {
        int pos = rng.randrange((int)pairs.size());
        String pair = pairs[pos];
        s1 += pair[0];
        s2 =  pair[1] + s2;
        ss1 += "(";
        ss2 += ")";
    }
    
    return Strings({s1 + "+" + s2, ss1 + "+" + ss2});
}

int
MotifTreeTopologyUnittest::test_creation() {
    //SS_Tree ss_t("(((+)))", "GAG+CUC");
    
    Strings h_seq_and_ss = random_helix();
    SS_Tree ss_t(h_seq_and_ss[1], h_seq_and_ss[0]);
    MotifTreeTopology mtt(ss_t);

    //std::cout << h_seq_and_ss[0] << std::endl;
    //std::cout << mtt.nodes().size() << std::endl;
    
    
    /*for(auto const & n : mtt.nodes()) {
        std::cout << n->data()->seq[0] << " " <<  n->data()->seq[1] << " " << n->data()->ends.size() << std::endl;
        std::cout << "end chains" << std::endl;
        for(auto const & end : n->data()->ends) {
            for(auto const & c : end.chains) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }
    }*/
    
    return 1;
}

int
MotifTreeTopologyUnittest::test_simple_build() {
    
    Strings h_seq_and_ss = random_helix();
    SS_Tree ss_t(h_seq_and_ss[1], h_seq_and_ss[0]);
    MotifTreeTopology mtt(ss_t);
    
    std::cout << h_seq_and_ss[0] << std::endl;
    MotifTree mt;
    
    
    for(auto const & n : mtt.nodes()) {
        //std::cout << n->data()->seq[0] << " " << n->data()->seq[1] << std::endl;
        //std::cout << name << std::endl;
        std::cout << "|" << n->data()->name << "|" << std::endl;
        MotifOP m  = LibraryManager::getInstance().get_motif(n->data()->name);
    }
    
    
    return 1;
}


int
MotifTreeTopologyUnittest::run() {
    if (test_creation() == 0)            { std::cout << "test_creation failed" << std::endl;  }
    if (test_simple_build() == 0)        { std::cout << "test_simple_build failed" << std::endl;  }
    return 0;
}


























