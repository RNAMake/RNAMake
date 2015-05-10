//
//  motif_tree_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_unittest.h"
#include "util/file_io.h"
#include "motif/motif_tree.h"
#include "resources/motif_library.h"


int
MotifTreeUnittest::test_creation() {
    MotifTree mt;
    return 1;
}


int
MotifTreeUnittest::test_add_motif() {
    MotifLibrary mlib (HELIX);
    MotifTree mt;
    MotifOP m = mlib.get_motif("HELIX.IDEAL");
    for (int i = 0; i < 20; i ++) { mt.add_motif(m); }
    
    MotifTreeNodeOP n = mt.nodes()[0];
    try {
        mt.add_motif(m, n);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(char const * e) { }
    
    try {
        mt.add_motif(m, nullptr, 10);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(char const * e) { }
    
    
    try {
        mt.add_motif(m, nullptr, -1, 10);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(char const * e) { }
    

    return 1;
}

int
MotifTreeUnittest::test_motif_tree_to_str() {
    String path = unittest_resource_dir() + "/motif_tree/test_motif_tree_to_str.dat";
    Strings lines = get_lines_from_file(path);
    ResidueTypeSet rts;
    MotifLibrary mlib ( TWOWAY );
    int i = 0;
    for (auto const & line : lines) {
        MotifTree mt ( line, rts );
        MotifTree mt2;
        int j = -1;
        for (auto const & n : mt.nodes()) {
            j++;
            if (j == 0) { continue; }
            MotifOP m = mlib.get_motif(n->motif()->name());
            MotifTreeNodeOP node = mt2.add_motif(m, nullptr, 0);
            if(node == nullptr) {
                std::cout << m->name() << std::endl;
            }
        }
        Point org = mt.nodes().back()->motif()->ends()[1]->d();
        Point np = mt2.nodes().back()->motif()->ends()[1]->d();
        float dist = org.distance(np);
        if (dist > 0.1) {
            return 0;
        }
        i++;
        
    }
    return 1;
}


int
MotifTreeUnittest::run() {
    if (test_creation() == 0)          { std::cout << "test_creation failed" << std::endl;  }
    if (test_add_motif() == 0)         { std::cout << "test_add_motif failed" << std::endl; }
    if (test_motif_tree_to_str() == 0) { std::cout << "test_motif_tree_to_str failed" << std::endl; }

    return 0;
}