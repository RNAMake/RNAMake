#ifndef __RNAMAKE_EXCEPTION_H__
#define __RNAMAKE_EXCEPTION_H__

#include <stdexcept>

#include <base/types.h>

namespace base {
//base project exception
struct RNAMakeException : public std::runtime_error {
    explicit RNAMakeException(String const& msg) : std::runtime_error(msg) {

    }
    explicit RNAMakeException(const char* msg) : std::runtime_error(msg) {

    }
};
//IO excpetion
struct RNAMakeIOException : public RNAMakeException {
    explicit RNAMakeIOException(String const& msg) : RNAMakeException(msg) {

    }

    explicit RNAMakeIOException(const char* msg) : RNAMakeException(msg) {

    }
};

//Implementation Exception
struct RNAMakeImplementationExcepetion : public RNAMakeException {
    explicit RNAMakeImplementationExcepetion(String const& msg) : RNAMakeException(msg) {

    }

    explicit RNAMakeImplementationExcepetion(const char* msg) : RNAMakeException(msg) {

    }
};

}

#endif // __RNAMAKE_EXCEPTION_H__
