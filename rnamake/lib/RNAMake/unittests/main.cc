//
//  main.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "atom_unittest.h"

int run_structure_unittests() {
    
    AtomUnittest atom_unittest;
    atom_unittest.run();
    
    return 1;
}



int main(int argc, const char * argv[]) {

    run_structure_unittests();


}