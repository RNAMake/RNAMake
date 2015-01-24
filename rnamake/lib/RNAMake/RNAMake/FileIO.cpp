//
//  FileIO.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/29/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "FileIO.h"

Strings
split_str_by_delimiter(
	String s,
	String delimiter) {
	
	std::string token;
	std::vector<std::string> tokens;
	size_t pos = 0;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		tokens.push_back(token);
		s.erase(0, pos + delimiter.length());
	}
	
	if(s.length() > 0) { tokens.push_back(s); }
	
	return tokens;
	
}

String
get_base_dir()  {
	char* base_path;
	base_path = std::getenv ("MotifAssembler");
	if (base_path==NULL) {
		std::cout << "cannot find environmental path MotifAssembler, please set it" << std::endl;
		exit(EXIT_FAILURE);
	}
	return String(base_path);
}