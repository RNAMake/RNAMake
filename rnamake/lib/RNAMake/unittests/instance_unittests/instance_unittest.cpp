//
//  instance_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/13/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "instance_unittest.hpp"
#include "instances/structure_instances.hpp"


namespace unittests {
namespace instances {

void
InstanceUnittests::test_residue() {
    auto r = ::instances::residue();
    std::cout << r->name() << std::endl;
}

int
InstanceUnittests::run() {
    test_residue();
    return 0;
}

int
InstanceUnittests::run_all() {
    return 0;
}

    
}
}