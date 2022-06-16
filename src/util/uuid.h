//
//  uuid.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_uuid_h
#define RNAMake_uuid_h

#include <fstream>
#include <iostream>

// RNAMake Headers
#include <base/types.hpp>

namespace util {

class Uuid {
public:
  Uuid();

  Uuid(Uuid const &uuid) : id_(uuid.id_) {}

  ~Uuid() {}

public:
  inline bool operator==(Uuid const &uuid) const { return id_ == uuid.id_; }

  inline bool operator!=(Uuid const &uuid) const { return id_ != uuid.id_; }

  inline bool operator<(Uuid const &uuid) const { return id_ < uuid.id_; }

  inline bool operator>(Uuid const &uuid) const { return id_ > uuid.id_; }

  friend std::ostream &operator<<(std::ostream &stream, Uuid const &uuid) {
    stream << uuid.id_;
    return stream;
  }

  uint64_t get_id() const { return id_; }

private:
  uint64_t id_;
};

typedef std::shared_ptr<Uuid> UuidOP;

struct UuidCompare {
  bool operator()(Uuid const &u1, Uuid const &u2) const {
    return u1.get_id() < u2.get_id();
  }
};

} // namespace util

#endif