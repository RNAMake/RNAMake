//
// Created by Joseph Yesselman on 2/14/17.
//

#ifndef RNAMAKE_SIMPLE_STRING_H
#define RNAMAKE_SIMPLE_STRING_H

#include "base/types.h"
#include "memory"
#include "sstream"
#include "iostream"

class SimpleString {
public:
    friend
    std::ostream & operator << (std::ostream &, SimpleString const &);

public:
    inline
    SimpleString() {
        chars_ = new char[1];
        size_ = 1;
    }

    inline
    SimpleString(String const & s) {
        chars_ = new char[s.length()];
        int i = 0;
        for(auto const & e : s) { chars_[i] = e; i++; }
        size_ = s.length();
    }


    inline
    SimpleString(SimpleString const & ss) {
        chars_ = ss.chars_;
        size_ = ss.size_;
    }

    ~SimpleString() {
        delete[] chars_;
    }

public:
    inline
    bool
    operator == (SimpleString const & other) const {
        for(int i = 0; i < size_; i++) {
            if(chars_[i] != other.chars_[i]) { return false; }
        }
        return true;
    }

    inline
    bool
    operator != (SimpleString const & other) const {
        for(int i = 0; i < size_; i++) {
            if(chars_[i] != other.chars_[i]) { return true; }
        }
        return false;
    }

public:
    String
    to_str() const {
        return String(chars_, size_);
    }

private:
    char * chars_;
    size_t size_;
};

std::ostream &
operator <<( std::ostream &, SimpleString const &);

typedef std::shared_ptr<SimpleString> SimpleStringOP;
typedef std::vector<SimpleStringOP>   SimpleStringOPs;

#endif //TEST_SIMPLE_STRING_H
