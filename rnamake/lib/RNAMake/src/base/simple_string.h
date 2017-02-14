//
// Created by Joseph Yesselman on 2/14/17.
//

#ifndef TEST_SIMPLE_STRING_H
#define TEST_SIMPLE_STRING_H

#include "types.h"
#include "memory"
#include "sstream"

class SimpleString {
public:
    friend
    std::ostream & operator << (std::ostream &, SimpleString const &);

public:
    inline
    SimpleString(String const & s) {
        chars_ = new char[s.length()];
        int i = 0;
        for(auto const & e : s) { chars_[i] = e; i++; }
        org_ = true;
    }

    inline
    SimpleString(SimpleString const & ss) {
        chars_ = ss.chars_;
        org_ = false;
    }

    inline
    ~SimpleString() {
        if(org_) { delete[] chars_; }
    }

private:
    char * chars_;
    bool org_;
};

std::ostream &
operator <<( std::ostream & stream, SimpleString const & ss) {
    auto string_stream = std::stringstream();
    string_stream << ss.chars_;
    stream << string_stream.str();
    return stream;
}


#endif //TEST_SIMPLE_STRING_H
