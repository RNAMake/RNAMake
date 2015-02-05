//
//  main.cpp
//  motif_tree.unittest
//
//  Created by Joseph Yesselman on 2/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "motif_type.h"
#include "motif_library.h"
#include "motif_tree.h"

int
test_creation() {
    MotifTree mt;
    mt.nodes()[0]->motif()->to_pdb("test.pdb");
    return 1;
}

int
test_add_motif() {
    MotifLibrary mlib (HELIX);
    MotifTree mt;
    MotifOP m = mlib.get_motif("HELIX.IDEAL");
    MotifOP m2 = mlib.get_motif("HELIX.IDEAL");
    mt.add_motif(m, NULL, 0, 0);
    MotifOP m3 = MotifOP ( new Motif(m->copy()));
    //std::cout << m->ends()[1]->r() << std::endl;
    //std::cout << m3->ends()[1]->r() << std::endl;
    mt.add_motif(m, NULL, 0, 1);
    mt.write_pdbs();
    
    
    return 1;
}




int main(int argc, const char * argv[]) {
    //if (test_creation() == 0)    { std::cout << "test_creation failed" << std::endl; }
    if (test_add_motif() == 0)   { std::cout << "test_add_motif failed" << std::endl; }
    return 0;
}
