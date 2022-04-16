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
  explicit Exception(String const &message) : std::runtime_error(message) {}
};

class ResourceException : public Exception {
public:
  explicit ResourceException(String const &message) : Exception(message) {}
};

class InputException : public Exception {
public:
  explicit InputException(String const &message) : Exception(message) {}
};

template<typename E>
inline void log_and_throw(String const & msg) {
  LOG_ERROR << msg;
  throw E(msg);
}

} // namespace base

#endif // RNAMAKE_SRC_BASE_EXCEPTION_HPP_
