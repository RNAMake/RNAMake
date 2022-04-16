
//
//  string.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//
#include <algorithm>
#include <functional>

// RNAMake Headers
#include <base/string.hpp>

namespace base::string {

// @brief splits string with a delimiter
// Known issue with escaped characters, will not work properly.
Strings split(String s, String const &delimiter) {
  String token;
  Strings tokens;
  size_t pos;
  while ((pos = s.find(delimiter)) != String::npos) {
    token = s.substr(0, pos);
    tokens.push_back(token);
    s.erase(0, pos + delimiter.length());
  }
  if (s.length() > 0) {
    tokens.push_back(s);
  }
  // to be consistent would be nice to always have at least 1 element
  if (tokens.empty()) {
    tokens.emplace_back();
  }
  return tokens;
}

// @brief joins a vector of strings with a delimiter
String join(Strings const &strs, String const &delimiter) {
  String return_s;
  int i = 0;
  for (auto const &s : strs) {
    return_s += s;
    if (i != strs.size() - 1) {
      return_s += delimiter;
    }
    i += 1;
  }
  return return_s;
}

// @brief Trims whitespace from the left end of the provided String
String left_trim(String s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
          }));
  return s;
}

// @brief Trims whitespace from the right end of the provided String
String right_trim(String s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); })
              .base(),
          s.end());
  return s;
}

// @brief Trims whitespace from both ends of the provided String
String trim(String s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
          }));
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); })
              .base(),
          s.end());
  return s;
}

// @brief adds quotes to around the string. string -> 'string'
String quoted(String const &s) { return String("'") + s + String("'"); }

}