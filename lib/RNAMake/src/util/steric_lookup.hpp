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

namespace util {

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
            math::Point const &);

    void
    add_points(
            math::Points const &);

    int
    clash(
            math::Point const &);

    int
    clash(
            math::Points const &);

    int
    better_clash(
            math::Point const &);

    int
    total_clash(
            math::Point const &);

    int
    total_clash(
            math::Points const &);

private:
    void
    _setup_additions();

private:
    std::map<double, int> bhash_;
    math::Points additions_, check_additions_;
    math::Point rounded_;
    math::Point p_;
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
            math::Point const &);

    void
    add_points(
            math::Points const &);

    bool
    clash(
            math::Point const &);

    bool
    clash(
            math::Points const &);

public:
    void
    to_pdb(
            String const &);

    int
    size() {
        return histo_.size(); }

private:
    void
    _setup_additions();

private:
    float grid_size_;
    float cutoff_;
    int radius_;
    math::Points additions_;
    math::ThreeDHistogram histo_;
    math::Point dummy_;
};

typedef std::shared_ptr<StericLookup> StericLookupOP;
typedef std::shared_ptr<StericLookupNew> StericLookupNewOP;

}

#endif /* steric_lookup_hpp */

























