//
// Created by Hassan Abdelsamad on 6/6/21.
//

#ifndef RNAMAKE_STRING_UTIL_H

#include <sstream>

std::string replace_char(std::string str, char ch1, char ch2) {
    for (int i = 0; i < str.length(); ++i) {
        if (str[i] == ch1)
            str[i] = ch2;
    }

    return str;
}

#define RNAMAKE_STRING_UTIL_H

#endif //RNAMAKE_STRING_UTIL_H
