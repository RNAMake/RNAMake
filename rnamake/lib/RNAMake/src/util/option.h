//
//  Option.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 11/9/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__Option__
#define __REDESIGNC__Option__

#include <stdio.h>
#include <vector>

//RNAMake Headers
#include "base/types.h"
#include "util/file_io.h"

struct Option {
	Option() {}
	Option(String nkey, String nvalue) { key = nkey; value = nvalue; }
	~Option() {}
	String key,value;
};

typedef std::vector<Option> Options;

Options
parse_command_into_options(
	int const,
	char const **,
    Strings allowed_options);

StringStringMap
parse_option_str(
	String);

StringStringMap
parse_file_options(
	String file_path);


#endif /* defined(__REDESIGNC__Option__) */
