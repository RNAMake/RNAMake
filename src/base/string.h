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

#include <cctype>

// RNAMake Headers
#include <base/types.h>

namespace base {

Strings split_str_by_delimiter(String, String);

Strings tokenize_line(String const &);

String join_by_delimiter(Strings const &, String const &);

String filename(String const &);

String base_dir(String const &);

DataType determine_string_data_type(String const &);

bool is_number(String const &);

String &ltrim(String &s);

String &rtrim(String &s);

String &trim(String &s);

String &replace_all(String &, String const &, String const &);
}  // namespace base

#endif /* defined(__RNAMake__string__) */