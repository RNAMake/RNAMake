
//
//  string.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__string__
#define __RNAMake__string__

#include <cstdio>

// RNAMake Headers
#include <base/types.hpp>

namespace base::string {

Strings split(String const &, const String &);

Strings tokenize(std::string const &, char);

String join(const Strings &, const String &);

String left_trim(String s);

String right_trim(String s);

String trim(String s);

String quoted(const String &);

template <typename... Args>
std::string format(const std::string &format, Args... args) {
  // Extra space for '\0'
  size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1;
  if (size <= 0) {
    throw std::runtime_error("Error during formatting.");
  }
  std::unique_ptr<char[]> buf(new char[size]);
  snprintf(buf.get(), size, format.c_str(), args...);
  return {buf.get(), buf.get() + size - 1}; // We don't want the '\0' inside
}

} // namespace base::string

#endif /* defined(__RNAMake__string__) */