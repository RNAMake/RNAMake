//
//  uuid_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>

#include "uuid_unittest.h"

int
UuidUnittest::test_compare() {
    Uuid u1, u2;
    
    if(u1 == u2) { return 0; }
    
    Uuid u3 = u1;
    
    if(!(u3 == u1)) { return 0; }
    
    return 1;
}
    
int
UuidUnittest::test_map() {
    
    std::map<Uuid, int, UuidCompare> uuid_int_map;
    
    Uuid u1;
    uuid_int_map[u1] = 1;

    for(int i = 0; i < 1000; i++) {
        Uuid u;
        uuid_int_map[u] = 1;
    }
    
    if(uuid_int_map.find(u1) == uuid_int_map.end()) {
        std::cout << "made it" << std::endl;
    }
    
    return 1;
    
}


int
UuidUnittest::run() {
    if (test_compare() == 0)    {  std::cout << "test_compare failed" << std::endl; }
    if (test_map() == 0)        {  std::cout << "test_map failed" << std::endl; }

    return 1;

}