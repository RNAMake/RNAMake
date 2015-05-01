//
//  motif_library_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_library_unittest.h"


int
MotifLibraryUnittest::test_load_all() {
    MotifLibrary mlib1(TWOWAY);
    mlib1.load_all();
    
    MotifLibrary mlib2(HELIX);
    mlib2.load_all();
    
    MotifLibrary mlib3(TCONTACT);
    mlib3.load_all();

    MotifLibrary mlib4(NWAY);
    mlib4.load_all();

    MotifLibrary mlib5(HAIRPIN);
    mlib5.load_all();
    return 1;
}

int
MotifLibraryUnittest::test_get_motif() {
    MotifLibrary mlib(HELIX);
    MotifOP m = mlib.get_motif("HELIX.IDEAL");

    try {
        MotifOP m2 = mlib.get_motif("HELIX.IDEALX");
        std::cout << "did not get exception" << std::endl;
        exit(0);
    } catch(String const & e) { }
    return 1;
}

int
MotifLibraryUnittest::test_contains_motif() {
    MotifLibrary mlib(HELIX);
    bool contains = mlib.contains_motif("HELIX.IDEAL");
    if(contains == false) { return 0; }
    bool contains2 = mlib.contains_motif("HELIX.IDEALX");
    if(contains2 == true) { return 0; }
    return 1;
}


int
MotifLibraryUnittest::run() {
    
    if (test_load_all() == 0)              { std::cout << "test_load_all failed" << std::endl; }
    if (test_get_motif() == 0)             { std::cout << "test_get_motif failed" << std::endl; }
    if (test_contains_motif() == 0)        { std::cout << "test_contains_motif failed" << std::endl; }
    
    return 0;
}