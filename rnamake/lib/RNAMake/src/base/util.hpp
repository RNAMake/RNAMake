//
//  util.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/22/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef base_util_h
#define base_util_h

#include <vector>
#include <map>
#include <algorithm>

#include "base/types.h"

template<typename T>
bool
element_in_vector(
        T const & element,
        std::vector<T> const & vec) {
    return std::find(vec.begin(), vec.end(), element) != vec.end();
}

bool
are_char_vectors_equal(
        Chars const & c1,
        Chars const & c2) {

    if(c1.size() != c2.size()) { return false; }
    for(int i = 0; i < c1.size(); i++) {
        if(c1[i] != c2[i]) { return false; }
    }
    return true;

}

String
char_vector_to_str(Chars const & c) {
    auto s = String();
    for(auto const & e : c) { s += e; }
    return s;
}

#endif /* base_util_h */
