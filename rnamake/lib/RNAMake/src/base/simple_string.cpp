//
// Created by Joseph Yesselman on 2/14/17.
//

#include "base/simple_string.h"

std::ostream &
operator <<( std::ostream & stream, SimpleString const & ss) {
    auto string_stream = std::stringstream();
    string_stream << ss.chars_;
    stream << string_stream.str();
    return stream;
}