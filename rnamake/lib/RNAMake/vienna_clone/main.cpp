//
//  main.cpp
//  vienna_clone
//
//  Created by Joseph Yesselman on 4/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "vienna_clone.h"

int main(int argc, const char * argv[]) {
    ViennaClone vc;
    vc.init_fold(100);
    vc.fold("AAAAAAAAAAUUUUUUUUUU");
    
    return 0;
}
