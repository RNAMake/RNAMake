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
#include "math/hashing.h"

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

class StericLookupNew {
public:
    StericLookupNew();

public:
    void
    add_point(
            Point const &);

    void
    add_points(
            Points const &);

    bool
    clash(
            Point const &);

    bool
    clash(
            Points const &);

public:
    void
    to_pdb(
            String const &);

private:
    void
    _setup_additions();

private:
    float grid_size_;
    float cutoff_;
    int radius_;
    Points additions_;
    ThreeDHistogram histo_;
    Point dummy_;
};

typedef std::shared_ptr<StericLookup> StericLookupOP;
typedef std::shared_ptr<StericLookupNew> StericLookupNewOP;

#endif /* steric_lookup_hpp */

























