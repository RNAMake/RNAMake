//
//  motif_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/9/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_connection__
#define __RNAMake__motif_connection__

#include <stdio.h>

#include "base/types.h"

class MotifConnection {
public:
    MotifConnection() {}
    
    inline
    MotifConnection(
        int i,
        int j,
        String const & name_i,
        String const & name_j):
    i_(i),
    j_(j),
    name_i_(name_i),
    name_j_(name_j)
    {}
    
    inline
    MotifConnection(
        MotifConnection const & mc):
    i_(mc.i_),
    j_(mc.j_),
    name_i_(mc.name_i_),
    name_j_(mc.name_j_)
    {}
    
    ~MotifConnection() {}
    
public:
    inline
    int
    i() { return i_; }
    
    inline
    int
    j() { return j_; }
    
    inline
    String const &
    name_i() { return name_i_; }
    
    inline
    String const &
    name_j() { return name_j_; }
    
private:
    int i_, j_;
    String name_i_, name_j_;
    
    
};

typedef std::shared_ptr<MotifConnection> MotifConnectionOP;
typedef std::vector<MotifConnectionOP>   MotifConnectionOPs;


#endif /* defined(__RNAMake__motif_connection__) */
