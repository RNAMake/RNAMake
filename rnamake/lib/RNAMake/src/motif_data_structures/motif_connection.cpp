//
//  motif_connection.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/9/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_connection.h"


String
MotifConnection::to_str() {
    String s = std::to_string(i_) + "," + std::to_string(j_) + "," + name_i_ + "," + name_j_;
    return s;
}