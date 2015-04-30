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

const char alphanum[] =
"0123456789"
"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
int stringLength = sizeof(alphanum) - 1;

Uuid::Uuid() {
    s_uuid_ = String();
    for ( int i = 0; i < 25; i++) {
        int pos = rand() % stringLength;
        s_uuid_ += alphanum[pos];
    }
    
}

std::ostream &
operator <<( std::ostream & stream, Uuid const & uuid) {
    stream << uuid.s_uuid();
    return stream;
}