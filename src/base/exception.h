#ifndef __RNAMAKE_EXCEPTION_H__
#define __RNAMAKE_EXCEPTION_H__

#include <stdexcept>

#include <base/types.h>

namespace base {
//base project exception
struct RNAMakeException : public std::runtime_error {
    RNAMakeException(String const& msg) : std::runtime_error(msg) {

    }
    RNAMakeException(const char* msg) : std::runtime_error(msg) {

    }
};
//IO excpetion
struct RNAMakeIOException : public RNAMakeException {
    RNAMakeIOException(String const& msg) : RNAMakeException(msg) {

    }

    RNAMakeIOException(const char* msg) : RNAMakeException(msg) {

    }
};

//Implementation Exception
struct RNAMakeImplementationExcepetion : public RNAMakeException {
    RNAMakeImplementationExcepetion(String const& msg) : RNAMakeException(msg) {

    }

    RNAMakeImplementationExcepetion(const char* msg) : RNAMakeException(msg) {

    }
};

}

#endif // __RNAMAKE_EXCEPTION_H__
