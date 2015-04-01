//
//  main.cpp
//  motif.unittest
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "settings.h"
#include "motif.h"
#include "residue_type_set.h"
#include "resource_manager.h"

Motif
get_test_motif() {
    String file = "test_str_to_motif.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Motif m(line, rts);
    return m;
}

int
test_str_to_motif() {
    String file = "test_str_to_motif.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Motif m(line, rts);
    return 1;
}

int
test_copy() {
    Motif m = get_test_motif();
    Motif mcopy = m.copy();
    Point p(10,20,20);
    m.move(p);
    m.to_pdb("original.pdb");
    mcopy.to_pdb("copy.pdb");
    return 1;
}

int
test_to_str() {
    Motif m = get_test_motif();
    ResidueTypeSet rts;
    Motif m2 = Motif(m.to_str(), rts);
    return 1;
}

int
test_update() {
    Motif m = get_test_motif();
    BasepairOP motif_end = m.basepairs()[0];
    std::cout << motif_end->res1()->get_atom("C1'")->coords() << std::endl;
    //std::cout << res1->atoms()[0]->coords() << std::endl;
    //std::cout << m.atoms()[0]->coords() << std::endl;
    //std::cout << m.residues()[0]->atoms()[0]->coords() << std::endl;
    m.move(Point(10,0,0));
    //ResidueOP res1 = m.get_residue(103, "A", "");
    //std::cout << res1->atoms()[0]->coords() << std::endl;
    std::cout << m.basepairs()[0]->res1()->get_atom("C1'")->coords() << std::endl;
    //motif_end = m.ends()[0];
    //std::cout << motif_end.res1()->get_atom("C1'")->coords() << std::endl;
    //std::cout << m.atoms()[0]->coords() << std::endl;
    //std::cout << m.residues()[0]->atoms()[0]->coords() << std::endl;
    
    return 1;
    
}

int
test_secondary_structure() {
    Motif m = get_test_motif();
    String org = "....((((((...(((.((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....)).)))..))))))(....((((..(((((((((.....)))))))).)))))......)";
    String new_ss = m.secondary_structure();
    //std::cout << new_ss << std::endl;
    for(int i = 0; i < org.size(); i++) {
        if(org[i] != new_ss[i]) { std::cout << i << " " << org[i] << " " << new_ss[i] << std::endl; }
    }
    
    return 1;
}

int
test_ref_motif() {
        
    for(int i = 0; i < 100000; i++) {
        Motif ref_m = ref_motif();
        std::cout << i << std::endl;
    }
    return 1;
}



int main(int argc, const char * argv[]) {
    //if (test_str_to_motif() == 0)         { std::cout << "test_str_to_motif failed" << std::endl; }
    //if (test_copy() == 0)                 { std::cout << "test_copy failed" << std::endl; }
    //if (test_to_str() == 0)               { std::cout << "test_to_str failed" << std::endl; }
    //if (test_update() == 0)               { std::cout << "test_update failed" << std::endl; }
    //if (test_secondary_structure() == 0)  { std::cout << "test_secondary_structure failed" << std::endl; }
    test_ref_motif();
    return 0;
}
