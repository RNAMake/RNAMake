//
//  motif_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "unittest.h"
#include "motif_unittest.h"
#include "util/file_io.h"


MotifUnittest::MotifUnittest() {
    String path = unittest_resource_dir() + "/motif/test_str_to_motif.dat";
    Strings lines = get_lines_from_file(path);
    
    ResidueTypeSet rts;
    m_ = Motif(lines[0], rts);
    
}

int
MotifUnittest::test_copy() {
    Motif mcopy = m_.copy();
    Point p(50,0,0);
    m_.move(p);
    
    Point org_center = center(m_.atoms());
    Point new_center = center(mcopy.atoms());
    float dist = (50) - org_center.distance(new_center);
    if(dist > 0.001) { return 0; }
    
    return 1;
}

int
MotifUnittest::test_to_str() {
    ResidueTypeSet rts;
    Motif m2 = Motif(m_.to_str(), rts);
    
    AtomOPs new_atoms = m2.atoms();
    AtomOPs org_atoms = m_.atoms();
    
    float dist;
    for(int i = 0; i < org_atoms.size(); i++) {
        if(org_atoms[i].get() == NULL) { continue; }
        dist = org_atoms[i]->coords().distance(new_atoms[i]->coords());
    }
    return 1;
}

int
MotifUnittest::test_secondary_structure() {
    String org = "....((((((...(((.((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....)).)))..))))))(....((((..(((((((((.....)))))))).)))))......)";
    String new_ss = m_.secondary_structure();
    for(int i = 0; i < org.size(); i++) {
        if(org[i] != new_ss[i]) { return 0; }
    }
    return 1;
}

int
MotifUnittest::run() {
    if (test_copy() == 0)                 { std::cout << "test_copy failed" << std::endl; }
    if (test_to_str() == 0)               { std::cout << "test_to_str failed" << std::endl; }
    if (test_secondary_structure() == 0)  { std::cout << "test_secondary_structure failed" << std::endl; }
    return 0;
}