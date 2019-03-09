//
//  steric_lookup.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef steric_lookup_hpp
#define steric_lookup_hpp

#include <map>
#include <stdio.h>

#include "math/xyz_vector.h"

class StericLookup {
public:
    StericLookup();

    StericLookup(
            float,
            float,
            int);
    
    ~StericLookup() {}

public:
    void
    add_point(
            Point const &);
    
    void
    add_points(
            Points const &);

    int
    clash(
            Point const &);
    
    int
    clash(
            Points const &);
    
    int
    better_clash(
            Point const &);

    int
    total_clash(
            Point const &);

    int
    total_clash(
            Points const &);
    
private:
    void
    _setup_additions();

private:
    std::map<double, int> bhash_;
    Points additions_, check_additions_;
    Point rounded_;
    Point p_;
    float grid_size_;
    float cutoff_;
    int radius_;
    double k_;
    
    
};

typedef std::shared_ptr<StericLookup> StericLookupOP;

#endif /* steric_lookup_hpp */
