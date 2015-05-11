//
//  string.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__string__
#define __RNAMake__string__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"

Strings
split_str_by_delimiter(
    String,
    String);

String
filename(
    String const &);

String
base_dir(
    String const &);


#endif /* defined(__RNAMake__string__) */
