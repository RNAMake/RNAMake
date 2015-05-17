//
//  chain_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__chain_unittest__
#define __RNAMake__chain_unittest__

#include <stdio.h>

#include "unittest.h"
#include "structure/chain.h"


class ChainUnittest : public Unittest {
public:
    ChainUnittest();
    
    ~ChainUnittest() {}
    

public:
    
    int
    test_to_str();
    
    int
    test_subchain();
    
    int
    test_to_pdb();
    
public:
    
    int
    run();
    
    void
    run_all();
    
private:
    ChainOP c_;
    
    
};

#endif /* defined(__RNAMake__chain_unittest__) */
