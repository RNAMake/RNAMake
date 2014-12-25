//
//  main.cpp
//  Atom.unittest
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "xyzVector.h"
#include "atom.h"

int
test_creation() {
    Atom atom("P", Point(0,1,2));
    std::cout << atom.coords() << std::endl;
    std::cout << atom.name() << std::endl;
    
    return 0;
}

int main(int argc, const char * argv[]) {
    test_creation();
    
    return 0;
}
