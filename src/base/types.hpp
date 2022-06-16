//
//  types.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_types_h
#define RNAMake_types_h

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

// renaming basic types to correspond with their use
using Index = int;
using Size = size_t;
using Flag = bool;
using String = std::string;
using Real = double;
using Path = std::filesystem::path;

// renaming vectors to reduce use of std::vector everywhere
using Indexes = std::vector<Index>;
using Chars = std::vector<char>;
using Reals = std::vector<Real>;
using Strings = std::vector<String>;

using StringOP = std::shared_ptr<String>;

typedef std::map<String, int> StringIntMap;
typedef std::map<String, float> StringFloatMap;
typedef std::map<String, String> StringStringMap;

// to keep track of data that may be many different types
// TODO unclear if I will use this?
enum class DataType {
  INT,
  FLOAT,
  STRING,
  BOOL,
  STRINGS,    // vector of string
  INTS,       // vector of ints
  FLOATS,     // vector of floats
  BOOLS,      // vector of bools
  XYZ_VECTOR, // {x,y,z} coords
  XYZ_MATRIX  // 3x3 matrix
};

typedef std::vector<DataType> DataTypes;

#endif
