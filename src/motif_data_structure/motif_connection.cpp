//
//  motif_connection.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/9/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structure/motif_connection.h"
#include "base/types.hpp"


namespace motif_data_structure {

String
MotifConnection::to_str() {
    String s = std::to_string(i_) + "," + std::to_string(j_) + "," + name_i_ + "," + name_j_;
    return s;
}

// motif connections constructors   ////////////////////////////////////////////////////////////////

MotifConnections::MotifConnections() :
        connections_(MotifConnectionOPs()) {}

MotifConnections::MotifConnections(
        MotifConnections const & mcs) :
        connections_(MotifConnectionOPs()) {

    for (auto const & c : mcs) {
        add_connection(c->i(), c->j(), c->name_i(), c->name_j());
    }

}

}