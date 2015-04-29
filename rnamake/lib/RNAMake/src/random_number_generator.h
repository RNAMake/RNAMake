//
//  random_number_generator.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_random_number_generator_h
#define RNAMake_random_number_generator_h
#include <random>

class RandomNumberGenerator {
public:
    RandomNumberGenerator() {
        srand(unsigned(time(NULL)));
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<float> dist(0,1);
        
        mt_ = mt;
        dist_ = dist;
    }
    
    
    inline
    float
    rand() { return dist_(mt_); }
    
    inline
    int
    randrange(int i) { return (int)(i*rand()); }
    
private:
    std::mt19937 mt_;
    std::uniform_real_distribution<float> dist_;
    
};


#endif
