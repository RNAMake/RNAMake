//
//  main.cpp
//  find_rings
//
//  Created by Joseph Yesselman on 3/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "cartesian_product.h"

int test_cartesian_product() {
    Ints int_test(3);
    int_test[0] = 0; int_test[1] = 1; int_test[2] = 2;
    std::vector<Ints> values(3);
    values[0] = int_test;
    values[1] = int_test;
    values[2] = int_test;
    
    CartesianProduct<int> test_product(values);
    while(! test_product.end() ) {
        Ints current = test_product.next();
        for(auto const & c : current) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }
    
    return 1;
}

int main(int argc, const char * argv[]) {
    // insert code here...
    
    std::cout << "Hello, World!\n";
    return 0;
}
