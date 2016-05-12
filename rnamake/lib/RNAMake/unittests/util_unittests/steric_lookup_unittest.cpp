//
//  steric_lookup_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "util/steric_lookup.hpp"
#include "steric_lookup_unittest.hpp"


Point const &
StericLookupUnittest::rand_point(int scale) {
    test_p_.x (rng_.rand() * rng_.randrange(scale));
    test_p_.y (rng_.rand() * rng_.randrange(scale));
    test_p_.z (rng_.rand() * rng_.randrange(scale));
    
    if(rng_.randrange(1000) < 500) { test_p_.x(-test_p_.x()); }
    if(rng_.randrange(1000) < 500) { test_p_.y(-test_p_.y()); }
    if(rng_.randrange(1000) < 500) { test_p_.z(-test_p_.z()); }

    return test_p_;
    
}

int
StericLookupUnittest::test_creation() {
    auto sl = StericLookup();
    
    return 0;
}

int
StericLookupUnittest::test_add_point() {
    auto sl = StericLookup();
    auto p = Point(0,0,0);
    sl.add_point(p);
    float grid_size_ = 0.5;
    auto additions_ = Points();
    
    auto add = Floats();
    for(int i = 1; i < 6; i++) {
        add.push_back(float(-i*grid_size_));
    }
    add.push_back(0);
    for(int i = 1; i < 6; i++) {
        add.push_back(float(i*grid_size_));
    }
    
    float dist = 0;
    Point p1;
    Point origin(0,0,0);
    for(auto const & x : add) {
        for(auto const & y : add) {
            for(auto const & z : add) {
                p1 = Point(x,y,z);
                dist = p1.distance(origin);
                if (dist > 3.5) {
                    additions_.push_back(p1);
                }
            }
        }
    }
    
    auto gp = Point();
    for(auto const & add : additions_) {
        gp = p + add;
        std::cout << sl.clash(gp) << " " << std::endl;
    
    }
    
    exit(0);
    auto rng = RandomNumberGenerator();

    auto test_p = Point(0,0,0);
    int clash = 0;
    for(int i = 0; i < 1000; i++) {
        clash = 0;
        test_p.x (rng.rand() * rng.randrange(10));
        test_p.y (rng.rand() * rng.randrange(10));
        test_p.z (rng.rand() * rng.randrange(10));
        
        dist = test_p.distance(p);
        if(dist < 2.5) {
            clash = 1;
        }
        
        std::cout << sl.clash(test_p) << " " << clash << " " << dist << std::endl;
        
    }
    
    return 0;
}

int
StericLookupUnittest::test_add_point_2() {
    auto sl = StericLookup();
    auto p = Point(0,0,0);
    sl.add_point(p);
    auto rng = RandomNumberGenerator();
    
    float dist = 0;
    auto test_p = Point(0,0,0);
    int clash = 0;
    for(int i = 0; i < 100000; i++) {
        clash = 0;
        test_p = rand_point();
        
        dist = test_p.distance(p);
        if(dist < 2.5) {
            clash = 1;
        }
        
        if(clash && ! sl.clash(test_p) ) {
            //std::cout << " made it" << std::endl;
        }
        
    }
    
    return 0;
}

int
StericLookupUnittest::test_add_point_3() {
    auto points = Points();
    for(int i = 0; i < 1000; i++) {
        points.push_back(rand_point(100));
    }
    
    auto sl = StericLookup();
    for (auto const & p : points) {
        sl.add_point(p);
    }
    
    float dist = 0;
    auto test_p = Point(0,0,0);
    int clash = 0;
    int false_count = 0;
    int miss_count = 0;
    int normal = 0;
    for(int i = 0; i < 100000; i++) {
        clash = 0;
        test_p = rand_point(100);
        
        for(auto const & p : points) {
            dist = test_p.distance(p);
            if(dist < 2.5) {
                clash = 1;
            }
        }
        
        if(clash && ! sl.clash(test_p) ) {
            miss_count += 1;
        }
        
        else if(sl.clash(test_p) == 1) {
            false_count += 1;
            
        }
        
        else {
            normal += 1;
        }
        
    }
    
    std::cout << false_count << " " << miss_count << " " << normal << std::endl;
    
    return 0;
}



int
StericLookupUnittest::run() {
    test_creation();
    //test_add_point();
    test_add_point_2();
    test_add_point_3();
    return 0;
}

int
StericLookupUnittest::run_all() {
    return 0;
}