//
//  uuid.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>
#include <iostream>

//RNAMake Headers
#include "util/uuid.h"
#include "util/random_number_generator.h"

namespace util {

const char alphanum[] =
"0123456789"
"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
int stringLength = sizeof(alphanum) - 1;

Uuid::Uuid() {
    s_uuid_ = String();
    RandomNumberGenerator rng;
    for ( int i = 0; i < 25; i++) {
        //int pos = rand() % stringLength;
        int pos = rng.randrange(stringLength);
        s_uuid_ += alphanum[pos];
    }
    
}

std::ostream &
operator <<( std::ostream & stream, Uuid const & uuid) {
    stream << uuid.s_uuid();
    return stream;
}

}