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

template<typename T>
bool
element_in_vector(T const & element, std::vector<T> const & vec) {
    return std::find(vec.begin(), vec.end(), element) != vec.end();
}


#endif /* base_util_h */
