//
//  types.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_types_h
#define RNAMake_types_h

#include <memory>
#include <string>
#include <vector>
#include <map>

typedef std::vector<int> Ints;
typedef std::vector<char> Chars;
typedef std::string String;
typedef std::vector<String> Strings;
typedef std::shared_ptr<String> StringOP;

typedef std::map<String, int> StringIntMap;
typedef std::map<String, float> StringFloatMap;
typedef std::map<String, String> StringStringMap;
#endif
