//
//  main.cpp
//  motif_tree.unittest
//
//  Created by Joseph Yesselman on 2/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <random>
#include "pose.h"
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
    for(int i = 0; i < 10; i++) {
        mt.add_motif(mlib.get_motif("HELIX.IDEAL"));
    }
    
    PoseOP pose = mt.to_pose();
    return 1;
}


int
test_merger2() {
    std::vector<MotifLibrary> mlibs;
    mlibs.push_back(ideal_helix_lib()); mlibs.push_back(unique_twoway_lib());
    MotifTree mt;
    int i = 0, pos = 0, count = 0, size = 100;
    int random = 0;
    //get crazy randomness
    srand(unsigned(time(NULL)));
    std::random_device rd;
    std::mt19937 mt19(rd());
    std::uniform_real_distribution<float> dist(0,1);
    while( i < size)  {
        if( i % 2 == 0) { pos = 0; }
        else            { pos = 1; }
        random = dist(mt19)*mlibs[pos].size();
        MotifOP m = mlibs[pos].motifs()[random];
        MotifTreeNodeOP node = mt.add_motif(m);
        if(node != NULL) { i++; }
        count++;
        if(count > 1000) { break; }
         
    }
    
    std::cout << mt.nodes().size() << std::endl;
    PoseOP pose = mt.to_pose();
    pose->to_pdb("merged.pdb");
    
    return 1;
}

int
test_memory_mangement() {
    MotifLibrary mlib (HELIX);
    for(int i = 0; i < 10000; i++) {
        MotifTree mt;
        mt.sterics(0);
        //MotifOP m = mlib.get_motif("HELIX.IDEAL");
        MotifOP m (new Motif(ref_motif()));
        for (int i = 0; i < 1; i ++) { mt.add_motif(m); }
        std::cout << i << std::endl;
        //break;
    }
    return 1;
}


int main(int argc, const char * argv[]) {
    //if (test_creation() == 0)          { std::cout << "test_creation failed" << std::endl;  }
    //if (test_add_motif() == 0)         { std::cout << "test_add_motif failed" << std::endl; }
    //if (test_motif_tree_to_str() == 0) { std::cout << "test_motif_tree_to_str failed" << std::endl; }
    //if (test_merger() == 0)            { std::cout << "test_merger failed" << std::endl;  }
    //if (test_merger2() == 0)            { std::cout << "test_merger2 failed" << std::endl;  }
    test_memory_mangement();
    return 0;
}













