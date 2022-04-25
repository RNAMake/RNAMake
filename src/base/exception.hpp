//
// Created by Joe Yesselman on 4/8/22.
//

#ifndef RNAMAKE_SRC_BASE_EXCEPTION_HPP_
#define RNAMAKE_SRC_BASE_EXCEPTION_HPP_

#include <stdexcept>

#include <base/types.hpp>
#include <base/log.hpp>

namespace base {

class Exception : public std::runtime_error {
public:
  explicit Exception(const String &message) : std::runtime_error(message) {}
};

class ResourceException : public Exception {
public:
  explicit ResourceException(const String &message) : Exception(message) {}
};

class InputException : public Exception {
public:
  explicit InputException(const String &message) : Exception(message) {}
};

class MathException : public Exception {
public:
  explicit MathException(const String &message) : Exception(message) {}
};

template<typename E>
inline void log_and_throw(const String & msg) {
  LOG_ERROR << msg;
  throw E(msg);
}

} // namespace base

#endif // RNAMAKE_SRC_BASE_EXCEPTION_HPP_
