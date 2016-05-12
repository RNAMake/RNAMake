//
//  steric_lookup.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <math.h>

#include "util/steric_lookup.hpp"


StericLookup::StericLookup():
bhash_(std::map<double, int>()),
grid_size_(0.5),
additions_(Points()){
    check_additions_ = Points();
    _setup_additions();
    
    
    
}

void
StericLookup::_setup_additions() {
    auto add = Floats();
    for(int i = 1; i < 6; i++) {
        add.push_back(float(-i*grid_size_));
    }
    add.push_back(0);
    for(int i = 1; i < 6; i++) {
        add.push_back(float(i*grid_size_));
    }
    
    float dist = 0;
    Point p;
    Point origin(0,0,0);
    for(auto const & x : add) {
        for(auto const & y : add) {
            for(auto const & z : add) {
                p = Point(x,y,z);
                dist = p.distance(origin);
                if (dist < 2.9 ) {
                    additions_.push_back(p);
                }
                if (dist <= 0.5) {
                    check_additions_.push_back(p);
                }
            }
        }
    }
}

void
StericLookup::add_point(
    Point const & p) {
    
    rounded_.x (round(p.x() / grid_size_)*grid_size_);
    rounded_.y (round(p.y() / grid_size_)*grid_size_);
    rounded_.z (round(p.z() / grid_size_)*grid_size_);
    
    auto gp = Point();
    double k;
    for(auto const & add : additions_) {
        gp = rounded_ + add;
        k = gp.x()*18397 + gp.y()*20483 + gp.z()*29303;
        bhash_[k] = 1;
    }

    
}

int
StericLookup::clash(
    Point const & p) {
    
    rounded_.x (round(p.x() / grid_size_)*grid_size_);
    rounded_.y (round(p.y() / grid_size_)*grid_size_);
    rounded_.z (round(p.z() / grid_size_)*grid_size_);
    
    double k = rounded_.x()*18397 + rounded_.y()*20483 + rounded_.z()*29303;
    
    if(bhash_.find(k) != bhash_.end()) { return 1; }
    else                              { return 0; }

    
}


int
StericLookup::better_clash(
    Point const & p) {
    
    rounded_.x (round(p.x() / grid_size_)*grid_size_);
    rounded_.y (round(p.y() / grid_size_)*grid_size_);
    rounded_.z (round(p.z() / grid_size_)*grid_size_);
    
    k_ = rounded_.x()*18397 + rounded_.y()*20483 + rounded_.z()*29303;
    
    if(bhash_.find(k_) == bhash_.end()) { return 0; }

    int i = 0;
    
    for(auto const & c_p : check_additions_) {
        p_ = rounded_ + c_p;
        k_ = p_.x()*18397 + p_.y()*20483 + p_.z()*29303;
        if(bhash_.find(k_) != bhash_.end()) { i++; }
    }
    

    if(i > 6) { return 1;}
    return 0;
    
    
    
    
}