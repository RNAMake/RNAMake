//
//  cartesian_product.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__cartesian_product__
#define __RNAMake__cartesian_product__

#include <stdio.h>
#include <iostream>

//RNAMake Headers
#include "base/types.h"

template< typename T >
class CartesianProduct {
public:
    typedef T Value;
    typedef std::vector<Value> Values;
    
public:
    CartesianProduct(std::vector<Values> const & values) {
        values_ = values;
        indices_ = Ints(values_.size());
        maxes_ = Ints(values_.size());
        current_ = Values(values_.size());
        end_ = 0;
        
        int i = 0;
        for (auto const & v : values_) {
            maxes_[i] = (int)v.size();
            i++;
        }
        
    }
    
    inline
    int const
    end() { return end_; }
    
    Values
    const &
    next() {
        int j = 0;
        for(auto const & v : values_) {
            current_ [j] = v [ indices_ [j] ];
            j++;
        }
        
        
        int i = (int)indices_.size()-1;
        while(i > -1) {
            indices_[i]++;
            
            if(indices_[i] == maxes_[i]) {
                if( i == 0) { end_ = 1; break; }
                indices_[i] = 0;
                i--;
                continue;
            }
            
            break;
        }
        return current_;
    }
    
private:
    Ints indices_;
    Ints maxes_;
    std::vector<Values> values_;
    Values current_;
    int end_;
    
};

#endif /* defined(__RNAMake__cartesian_product__) */
