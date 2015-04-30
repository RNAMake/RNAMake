//
//  main.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util_unittests/uuid_unittest.h"

#include "structure_unittests/atom_unittest.h"
#include "structure_unittests/residue_type_set_unittest.h"
#include "structure_unittests/residue_unittest.h"
#include "structure_unittests/chain_unittest.h"


int main(int argc, const char * argv[]) {

    ChainUnittest test;
    test.run();
    

}