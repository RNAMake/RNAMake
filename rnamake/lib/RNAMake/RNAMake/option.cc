//
//  Option.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 11/9/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "Option.h"
#include <fstream>

Options
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

StringStringMap
parse_option_str(
	String opt_string) {
	StringStringMap options;
	
	Strings spl = split_str_by_delimiter(opt_string,",");
	for(auto opt_and_value : spl) {
		Strings opt_strs = split_str_by_delimiter(opt_and_value,"=");
		if(opt_strs.size() == 2) {
			options[opt_strs[0]] = opt_strs[1];
		}
	}
	
	return options;
}


StringStringMap
parse_file_options(
	String file_path) {
	StringStringMap options;
	String line;
	std::ifstream in;
	in.open(file_path);
	while(in.good()) {
		getline(in,line);
		StringStringMap new_options = parse_option_str(line);
		for(auto const & kv : new_options) {
			options[kv.first] = kv.second;
		}
		
	}
	in.close();
	return options;
}




