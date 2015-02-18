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
#include "residue_type_set.h"

Strings
get_lines_from_file(String const fname) {
    String line;
    Strings lines;
    std::ifstream input;
    input.open(fname);
    while ( input.good() ) {
        getline(input, line);
        if( line.length() < 10 ) { break; }
        lines.push_back(line);
        
    }
    return lines;

}

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
    for (int i = 0; i < 20; i ++) { mt.add_motif(m); }
    //mt.write_pdbs();
    return 1;
}

int
test_motif_tree_to_str() {
    Strings lines = get_lines_from_file("test_motif_tree_to_str.dat");
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
            MotifTreeNodeOP node = mt2.add_motif(m, NULL, 0);
            if(node == NULL) {
                std::cout << m->name() << std::endl;
            }
        }
        Point org = mt.nodes().back()->motif()->ends()[1]->d();
        Point np = mt2.nodes().back()->motif()->ends()[1]->d();
        float dist = org.distance(np);
        if (dist > 0.1) {
            std::cout << i << " " << dist << " " << mt.nodes().size() << " " << mt2.nodes().size() << std::endl;
        }
        i++;
        
    }
    return 1;
}

int
test_merger() {
    MotifLibrary mlib(HELIX);
    MotifTree mt;
    mt.add_motif(mlib.get_motif("HELIX.IDEAL"));
    mt.add_motif(mlib.get_motif("HELIX.IDEAL"));
    mt.to_pose();
    return 1;
}

int main(int argc, const char * argv[]) {
    //if (test_creation() == 0)          { std::cout << "test_creation failed" << std::endl;  }
    //if (test_add_motif() == 0)         { std::cout << "test_add_motif failed" << std::endl; }
    //if (test_motif_tree_to_str() == 0) { std::cout << "test_motif_tree_to_str failed" << std::endl; }
    if (test_merger() == 0)            { std::cout << "test_merger failed" << std::endl;  }
    return 0;
}













