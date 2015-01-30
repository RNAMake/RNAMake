//
//  main.cpp
//  motif.unittest
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "motif.h"
#include "residue_type_set.h"

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


int main(int argc, const char * argv[]) {
    // insert code here...
    test_str_to_motif();
    std::cout << "Hello, World!\n";
    return 0;
}
