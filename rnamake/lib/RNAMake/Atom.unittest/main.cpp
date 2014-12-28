//
//  main.cpp
//  Atom.unittest
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "gtest/gtest.h"
#include <iostream>
#include "xyzVector.h"
#include "atom.h"
#include "numerical.h"


int
test_creation() {
    Atom atom("P", Point(0, 1, 2));
    Point p(0, 1, 2);
    std::cout << atom.to_pdb_str(1) << std::endl;
    if (!are_xyzVector_equal(atom.coords(), p)) {
        return 0;
    }
    return 1;
}

int main(int argc, const char * argv[]) {
    if (test_creation() == 0) {
        std::cout << "test_creation failed" << std::endl;
        exit(0);
    }
    
    

}
