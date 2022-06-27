//
//  uuid.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_uuid_h
#define RNAMake_uuid_h

#include <fstream>
#include <iostream>

// RNAMake Headers
#include <base/types.hpp>
#include <external/sole/sole.hpp>

namespace util {
typedef sole::uuid Uuid;
Uuid generate_uuid() { return sole::uuid0(); }

Uuid uuid_from_str(const String &str) { return sole::rebuild(str); }

} // namespace util

#endif