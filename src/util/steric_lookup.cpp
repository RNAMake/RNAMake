//
//  steric_lookup.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <math.h>

#include "util/steric_lookup.hpp"
#include <math/hashing.h>
#include <math/matrix_3x3.hpp>
#include <math/vector_3.hpp>

namespace util {

// constructors
/*
/// @brief - creates a map of a certain grid size
// this is a constructor
StericLookup::StericLookup()
    : _bhash(std::map<double, int>()), _grid_size(0.5), _cutoff(2.65),
      _radius(6), _additions(math::Vector3s()) {
  _check_additions = math::Vector3s();
  _setup_additions();
}

/// @brief - sets up a grid of unit size grid_size, total grid size cutoff, and
/// buffer radius radius
/// note: this is a constructor
StericLookup::StericLookup(float grid_size, float cutoff, int radius)
    : _bhash(std::map<double, int>()), _grid_size(grid_size), _cutoff(cutoff),
      _radius(radius), _additions(math::Vector3s()) {
  _check_additions = math::Vector3s();
  _setup_additions();
}

/// @brief -
void StericLookup::_setup_additions() {
  auto add = Reals(); // initializes a vector of real numbers
  for (int i = 1; i < _radius; i++) {
    add.push_back(float(-i * _grid_size));
  }
  add.push_back(0);
  for (int i = 1; i < _radius; i++) {
    add.push_back(float(i * _grid_size));
  }

  float dist = 0;
  math::Vector3 p;
  math::Vector3 origin(0, 0, 0);
  for (auto const &x : add) {
    for (auto const &y : add) {
      for (auto const &z : add) {
        p = math::Vector3(x, y, z);
        dist = p.distance(origin);
        if (dist < _cutoff) {
          _additions.push_back(p);
        }
      }
    }
  }
}

/// @brief - adds a point in space for an atom
void StericLookup::add_point(math::Vector3 const &p) {
  _rounded.set_x(round(p.get_x() / _grid_size) * _grid_size);
  _rounded.set_y(round(p.get_y() / _grid_size) * _grid_size);
  _rounded.set_z(round(p.get_z() / _grid_size) * _grid_size);

  auto gp = math::Vector3();
  double k;

  for (auto const &add : _additions) {
    gp = _rounded + add;
    k = double(gp.get_x() * 0.001) + (double)(gp.get_y()) * 0.000000001 +
        double(gp.get_z() * 100000);

    if (_bhash.find(k) != _bhash.end()) {
      _bhash[k] += 1;
    } else {
      _bhash[k] = 1;
    }
  }
}

/// @brief - adds points for atoms from an array of (position) vectors
void StericLookup::add_points(math::Vector3s const &points) {
  for (auto const &p : points) {
    add_point(p);
  }
}

/// @brief - checks if a given point clashes with/has an overlapping radius with
/// other points
bool StericLookup::clash(math::Vector3 const &p) {
  _rounded.set_x(round(p.get_x() / _grid_size) * _grid_size);
  _rounded.set_y(round(p.get_y() / _grid_size) * _grid_size);
  _rounded.set_z(round(p.get_z() / _grid_size) * _grid_size);

  double k = double(_rounded.get_x() * 0.001) +
             (double)(_rounded.get_y()) * 0.000000001 +
             double(_rounded.get_z() * 100000);

  if (_bhash.find(k) != _bhash.end()) {
    return true;
  } else {
    return false;
  }
}

/// @brief - checks if a list of given points clashes with/has an overlapping
/// radius with other points
// TODO write a unittest
bool StericLookup::clash(math::Vector3s const &points) {
  bool is_clash = false;
  for (auto const &p : points) {
    is_clash = clash(p);
    if (is_clash) {
      return is_clash;
    }
  }
  return 0;
}
/// @brief -
int StericLookup::better_clash(math::Vector3 const &p) {
  _rounded.set_x(round(p.get_x() / _grid_size) * _grid_size);
  _rounded.set_y(round(p.get_y() / _grid_size) * _grid_size);
  _rounded.set_z(round(p.get_z() / _grid_size) * _grid_size);

  _k = _rounded.get_x() * 18397 + _rounded.get_y() * 20483 +
       _rounded.get_z() * 29303;

  if (_bhash.find(_k) == _bhash.end()) {
    return 0;
  }

  int i = 0;

  for (auto const &c_p : _check_additions) {
    _p = _rounded + c_p;
    _k = _p.get_x() * 18397 + _p.get_y() * 20483 + _p.get_z() * 29303;
    if (_bhash.find(_k) != _bhash.end()) {
      i++;
    }
  }

  if (i > 6) {
    return 1;
  }
  return 0;
}

/// @brief - counts the number of clashes in a lookup consisting of a single
/// vector
int StericLookup::total_clash(math::Vector3 const &p) {
  _rounded.set_x(round(p.get_x() / _grid_size) * _grid_size);
  _rounded.set_y(round(p.get_y() / _grid_size) * _grid_size);
  _rounded.set_z(round(p.get_z() / _grid_size) * _grid_size);

  double k = double(_rounded.get_x() * 0.001) +
             (double)(_rounded.get_y()) * 0.000000001 +
             double(_rounded.get_z() * 100000);

  if (_bhash.find(k) != _bhash.end()) {
    return _bhash[k];
  } else {
    return 0;
  }
}

/// @brief - counts the number of clashes in a lookup consititing of an array of
/// vectors
int StericLookup::total_clash(math::Vector3s const &points) {
  int clash_count = 0;
  for (auto const &p : points) {
    clash_count += clash(p);
  }

  return clash_count;
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// StericLookupNew - TODO make this the old one now
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StericLookupNew::StericLookupNew() {
  auto bb = math::BoundingBox(math::Vector3(-200, -200, -200),
                              math::Vector3(100, 100, 100));
  // auto bb = math::BoundingBox(math::Point(-10, -10, -10), math::Point(10, 10,
  // 10));
  auto bin_widths = math::Real3{0.25, 0.25, 0.25};
  _histo = math::ThreeDHistogram(bb, bin_widths);
  _grid_size = 0.25;
  _cutoff = 2.70;
  _radius = 12;
  _additions = math::Vector3s();
  _setup_additions();
}

void StericLookupNew::_setup_additions() {
  auto add = Reals();
  for (int i = 1; i < _radius; i++) {
    add.push_back(float(-i * _grid_size));
  }
  add.push_back(0);
  for (int i = 1; i < _radius; i++) {
    add.push_back(float(i * _grid_size));
  }

  float dist = 0;
  math::Vector3 p;
  math::Vector3 origin(0, 0, 0);
  for (auto const &x : add) {
    for (auto const &y : add) {
      for (auto const &z : add) {
        p = math::Vector3(x, y, z);
        dist = p.distance(origin);
        if (dist < _cutoff) {
          _additions.push_back(p);
        }
      }
    }
  }
}

/// @brief - adds a point to the steric lookup
void StericLookupNew::add_point(math::Vector3 const &p) {
  for (auto const &add : _additions) {
    if (p.get_x() < -200 || p.get_y() < -200 || p.get_z() < -200) {
      String msg = "Point is outside the boundary! Bounds are from -200 to 100";
      base::log_and_throw<base::MathException>(msg);
    }
    if (p.get_x() > 100 || p.get_y() > 100 || p.get_z() > 100) {
      String msg = "Point is outside the boundary! Bounds are from -200 to 100";
      base::log_and_throw<base::MathException>(msg);
    }
    _dummy = p + add;
    _histo.add(_dummy);
  }
}

/// @brief - adds a set of points to the steric lookup
void StericLookupNew::add_points(math::Vector3s const &points) {
  for (auto const &p : points) {
    add_point(p);
  }
}

/// @brief - checks if a point clashes with the lookup
bool StericLookupNew::clash(math::Vector3 const &p) {
  return _histo.contains(p);
}

/// @brief - checks if any points in a set clash with the lookup
bool StericLookupNew::clash(math::Vector3s const &points) {
  for (auto const &p : points) {
    bool is_clash = clash(p);
    if (is_clash == true) {
      return is_clash;
    }
  }
  return false;
}

void StericLookupNew::to_pdb(String const &pdb_name) {
  _histo.write_histo_to_pdb(pdb_name);
}
} // namespace util
