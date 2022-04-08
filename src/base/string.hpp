
//
//  string.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__string__
#define __RNAMake__string__

#include <cstdio>

// RNAMake Headers
#include <base/types.hpp>

namespace base::string {

Strings split(String, String const &);

String join(Strings const &, String const &);

String left_trim(String s);

String right_trim(String s);

String trim(String s);

String quoted(String const &);

}

#endif /* defined(__RNAMake__string__) */