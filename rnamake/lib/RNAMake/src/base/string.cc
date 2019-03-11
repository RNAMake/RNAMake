//
//  string.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//
#include <algorithm>
#include <functional>

//RNAMake Headers
#include "base/string.h"

namespace base {

Strings
split_str_by_delimiter(
        String s,
        String delimiter) {

    String token;
    std::vector<String> tokens;
    size_t pos = 0;
    while ((pos = s.find(delimiter)) != String::npos) {
        token = s.substr(0, pos);
        tokens.push_back(token);
        s.erase(0, pos + delimiter.length());
    }

    if (s.length() > 0) { tokens.push_back(s); }

    return tokens;

}

String
join_by_delimiter(
        Strings const & strs,
        String const & delimiter) {

    String return_s;
    int i = 0;
    for (auto const & s : strs) {
        return_s += s;
        if (i != strs.size() - 1) { return_s += delimiter; }
    }

    return return_s;


}


String
filename(
        String const & path) {
    Strings path_spl = base::split_str_by_delimiter(path, "/");
    return path_spl.back();
}


String
base_dir(
        String const & path) {
    Strings path_spl = base::split_str_by_delimiter(path, "/");
    String base_path;
    for (int i = 0; i < path_spl.size() - 1; i++) {
        base_path += path_spl[i] + "/";
    }
    return base_path;

}


bool
is_number(
        String const & s) {

    for (auto const & c : s) {
        if (!std::isdigit(c)) { return false; }
    }
    return true;

    //return !s.empty() && std::find_if(s.begin(),
    //                                  s.end(), [](char c) { return !std::isdigit(c); }) == s.end();

}

/**
 * @brief Left Trim
 *
 * Trims whitespace from the left end of the provided String
 *
 * @param[out] s The String to trim
 *
 * @return The modified String&
 */
String & ltrim(String & s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::ptr_fun<int, int>(std::isgraph)));
    return s;
}

/**
 * @brief Right Trim
 *
 * Trims whitespace from the right end of the provided String
 *
 * @param[out] s The String to trim
 *
 * @return The modified String&
 */
String & rtrim(String & s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::ptr_fun<int, int>(std::isgraph)).base(), s.end());
    return s;
}

/**
 * @brief Trim
 *
 * Trims whitespace from both ends of the provided String
 *
 * @param[out] s The String to trim
 *
 * @return The modified String&
 */
String & trim(String & s) {
    return ltrim(rtrim(s));
}

}




