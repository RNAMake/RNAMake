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
join_by_delimiter(
    Strings const &,
    String const &);

String
filename(
    String const &);

String
base_dir(
    String const &);

bool
is_number(
    String const &);

String & ltrim(
    String & s);

String& rtrim(
    String & s);

String& trim(
    String & s);

#endif /* defined(__RNAMake__string__) */
