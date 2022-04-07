
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

// RNAMake Headers
#include <base/types.hpp>

namespace base::string {

Strings split(String, String const &);

String join_by_delimiter(Strings const &, String const &);

String filename(String const &);

String base_dir(String const &);

bool is_number(String const &);

String &ltrim(String &s);

String &rtrim(String &s);

String &trim(String &s);

bool is_char_in_string(char, String const &);

String quoted_string(String const &);

String string_map_to_string(StringStringMap const &);

String &replace_all(String &, String const &, String const &);

} // namespace base

#endif /* defined(__RNAMake__string__) */