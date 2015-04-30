//
//  chain.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_chain_fwd_h
#define RNAMake_chain_fwd_h
#include <vector>
#include <memory>

class Chain;

typedef std::shared_ptr<Chain> ChainOP;
typedef std::vector<ChainOP> ChainOPs;

#endif
