//
// Created by Joe Yesselman on 6/14/22.
//

#ifndef RNAMAKE_SRC_UTIL_SQLITE_EXCEPTION_HPP_
#define RNAMAKE_SRC_UTIL_SQLITE_EXCEPTION_HPP_

#include <base/exception.hpp>

namespace util {

class SqliteException : public std::runtime_error {
public:
  explicit SqliteException(const String &message)
      : std::runtime_error(message) {}
};

} // namespace util
#endif // RNAMAKE_SRC_UTIL_SQLITE_EXCEPTION_HPP_
