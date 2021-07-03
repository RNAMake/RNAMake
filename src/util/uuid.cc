
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

namespace util {

    Uuid::Uuid() {
        int length = 10;
        id_ = rand() % length;
        uint64_t digit = 10;
        for (int i = 1; i < 15; i++) {
            int pos = rand() % length;
            id_ += (digit * i) * pos;
            digit *= 10;
        }
    }
}

/*std::ostream &
operator <<( std::ostream & stream, Uuid const & uuid) {
    stream << uuid.s_uuid();
    return stream;
}*/