//
//  Option.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 11/9/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include <fstream>

//RNAMake Headers
#include "base/string.h"
#include "base/option.h"

/*Options
parse_command_into_options(
	int const argc,
	char const ** argv,
    Strings allowed_options) {
	
	Options opts;
	int last_arg = argc-1;
	for(int i = 1; i < last_arg; i++) {
		if(argv[i][0] == '-') {
			String key(argv[i]);
			key = key.substr(1);
			
            //if (i < last_arg-1)                         { opts.push_back(Option(key,"")); continue; }
			if (i < last_arg && argv[i+1][0] != '-')	{ opts.push_back(Option(key,argv[i+1])); }
			else									    { opts.push_back(Option(key,"")); }
			
		}
	}
	return opts;
}
}*/




