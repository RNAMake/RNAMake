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

#include <math/hashing.h>
#include <math/vector_3.hpp>

namespace util {

class StericLookupNew {
public:
  StericLookupNew();

  StericLookupNew(float grid_size, float cutoff, int radius);

  ~StericLookupNew() {}

public:
  void add_point(math::Vector3 const &);

  void add_points(math::Vector3s const &);

  bool clash(math::Vector3 const &);

  bool clash(math::Vector3s const &);

//  int total_clash(math::Vector3 const &);

  int total_clash(math::Vector3s const &points);

public:
  void to_pdb(String const &);

  int size() { return _histo.size(); }

private:
  void _setup_additions();

private:
  float _grid_size;
  float _cutoff;
  int _radius;
  math::Vector3s _additions;
  math::ThreeDHistogram _histo;
  math::Vector3 _dummy;
};

// typedef std::shared_ptr<StericLookup> StericLookupOP;
typedef std::shared_ptr<StericLookupNew> StericLookupNewOP;

} // namespace util

#endif /* steric_lookup_hpp */
