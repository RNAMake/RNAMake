//
// Created by Joseph Yesselman on 10/22/17.
//

#ifndef RNAMAKE_NEW_ASSERTIONS_H
#define RNAMAKE_NEW_ASSERTIONS_H

#include <exception>
#include <memory>

#include <base/types.h>

template <typename exception>
void
expects(
    bool condition,
    String const & message) {

    if(!condition) { throw exception(message); }
}

template <typename exception>
void
ensures(
        bool condition,
        String const & message) {

    if (!condition) { throw exception(message); }
}



#endif //RNAMAKE_NEW_ASSERTIONS_H
