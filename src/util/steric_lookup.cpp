//
//  steric_lookup.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <math.h>

#include "util/steric_lookup.hpp"
#include <iomanip>      // std::setprecision

namespace util {

StericLookup::StericLookup() :
        bhash_(std::map<double, int>()),
        grid_size_(0.5),
        cutoff_(2.65),
        radius_(6),
        additions_(math::Points()) {
    check_additions_ = math::Points();
    _setup_additions();
}

StericLookup::StericLookup(
        float grid_size,
        float cutoff,
        int radius) :
        bhash_(std::map<double, int>()),
        grid_size_(grid_size),
        cutoff_(cutoff),
        radius_(6),
        additions_(math::Points()) {
    check_additions_ = math::Points();
    _setup_additions();
}

void
StericLookup::_setup_additions() {
    auto add = Floats();
    for (int i = 1; i < radius_; i++) {
        add.push_back(float(-i * grid_size_));
    }
    add.push_back(0);
    for (int i = 1; i < radius_; i++) {
        add.push_back(float(i * grid_size_));
    }

    float dist = 0;
    math::Point p;
    math::Point origin(0, 0, 0);
    for (auto const & x : add) {
        for (auto const & y : add) {
            for (auto const & z : add) {
                p = math::Point(x, y, z);
                dist = p.distance(origin);
                if (dist < cutoff_) {
                    additions_.push_back(p);
                }
            }
        }
    }
}

void
StericLookup::add_point(
        math::Point const & p) {

    rounded_.set_x(round(p.get_x() / grid_size_) * grid_size_);
    rounded_.set_y(round(p.get_y() / grid_size_) * grid_size_);
    rounded_.set_z(round(p.get_z() / grid_size_) * grid_size_);

    auto gp = math::Point();
    double k;

    for (auto const & add : additions_) {
        gp = rounded_ + add;
        k = double(gp.get_x() * 0.001) + (double) (gp.get_y()) * 0.000000001 + double(gp.get_z() * 100000);

        if (bhash_.find(k) != bhash_.end()) { bhash_[k] += 1; }
        else { bhash_[k] = 1; }
    }


}


void
StericLookup::add_points(
        math::Points const & points) {

    for (auto const & p : points) {
        add_point(p);
    }

}

int
StericLookup::clash(
        math::Point const & p) {

    rounded_.set_x(round(p.get_x() / grid_size_) * grid_size_);
    rounded_.set_y(round(p.get_y() / grid_size_) * grid_size_);
    rounded_.set_z(round(p.get_z() / grid_size_) * grid_size_);

    double k = double(rounded_.get_x() * 0.001) + (double) (rounded_.get_y()) * 0.000000001 + double(rounded_.get_z() * 100000);


    if (bhash_.find(k) != bhash_.end()) { return 1; }
    else { return 0; }

}

int
StericLookup::clash(
        math::Points const & points) {

    int is_clash = 0;
    for (auto const & p : points) {
        is_clash = clash(p);
        if (is_clash) { return is_clash; }
    }

    return 0;

}


int
StericLookup::better_clash(
        math::Point const & p) {

    rounded_.set_x(round(p.get_x() / grid_size_) * grid_size_);
    rounded_.set_y(round(p.get_y() / grid_size_) * grid_size_);
    rounded_.set_z(round(p.get_z() / grid_size_) * grid_size_);

    k_ = rounded_.get_x() * 18397 + rounded_.get_y() * 20483 + rounded_.get_z() * 29303;

    if (bhash_.find(k_) == bhash_.end()) { return 0; }

    int i = 0;

    for (auto const & c_p : check_additions_) {
        p_ = rounded_ + c_p;
        k_ = p_.get_x() * 18397 + p_.get_y() * 20483 + p_.get_z() * 29303;
        if (bhash_.find(k_) != bhash_.end()) { i++; }
    }


    if (i > 6) { return 1; }
    return 0;

}


int
StericLookup::total_clash(
        math::Point const & p) {
    rounded_.set_x(round(p.get_x() / grid_size_) * grid_size_);
    rounded_.set_y(round(p.get_y() / grid_size_) * grid_size_);
    rounded_.set_z(round(p.get_z() / grid_size_) * grid_size_);

    double k = double(rounded_.get_x() * 0.001) + (double) (rounded_.get_y()) * 0.000000001 + double(rounded_.get_z() * 100000);


    if (bhash_.find(k) != bhash_.end()) { return bhash_[k]; }
    else { return 0; }


}

int
StericLookup::total_clash(
        math::Points const & points) {
    int clash_count = 0;
    for (auto const & p : points) {
        clash_count += clash(p);
    }

    return clash_count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// StericLookupNew
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StericLookupNew::StericLookupNew() {
    auto bb = math::BoundingBox(math::Point(-200, -200, -200), math::Point(100, 100, 100));
    //auto bb = math::BoundingBox(math::Point(-10, -10, -10), math::Point(10, 10, 10));
    auto bin_widths = math::Real3{0.25, 0.25, 0.25};
    histo_ = math::ThreeDHistogram(bb, bin_widths);
    grid_size_ = 0.25;
    cutoff_ = 2.70;
    radius_ = 12;
    additions_ = math::Points();
    _setup_additions();
}

void
StericLookupNew::_setup_additions() {
    auto add = Floats();
    for (int i = 1; i < radius_; i++) {
        add.push_back(float(-i * grid_size_));
    }
    add.push_back(0);
    for (int i = 1; i < radius_; i++) {
        add.push_back(float(i * grid_size_));
    }

    float dist = 0;
    math::Point p;
    math::Point origin(0, 0, 0);
    for (auto const & x : add) {
        for (auto const & y : add) {
            for (auto const & z : add) {
                p = math::Point(x, y, z);
                dist = p.distance(origin);
                if (dist < cutoff_) {
                    additions_.push_back(p);
                }
            }
        }
    }
}

void
StericLookupNew::add_point(
        math::Point const & p) {
    for (auto const & add : additions_) {
        dummy_ = p + add;
        histo_.add(dummy_);
    }
}

void
StericLookupNew::add_points(
        math::Points const & points) {
    for (auto const & p : points) { add_point(p); }
}

bool
StericLookupNew::clash(
        math::Point const & p) {
    return histo_.contains(p);
}

bool
StericLookupNew::clash(
        math::Points const & points) {

    bool is_clash = 0;
    for (auto const & p : points) {
        is_clash = clash(p);
        if (is_clash) { return is_clash; }
    }

    return false;
}

void
StericLookupNew::to_pdb(
        String const & pdb_name) {
    histo_.write_histo_to_pdb(pdb_name);
}

}
























