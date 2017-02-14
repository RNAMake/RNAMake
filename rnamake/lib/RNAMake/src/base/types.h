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
#include <array>

using Ints = std::vector<int>;
using Chars = std::vector<char>;
using Shorts = std::vector<short>;
using Floats = std::vector<float>;

using String = std::string;
using StringOP = std::shared_ptr<String>;

using Strings = std::vector<String>;

typedef std::map<String, int> StringIntMap;
typedef std::map<String, float> StringFloatMap;
typedef std::map<String, String> StringStringMap;

template<typename T>
using   VectorUP = std::unique_ptr<std::vector<T>>;

#endif
