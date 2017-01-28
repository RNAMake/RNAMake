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

//RNAMake Headers
#include "base/types.h"


class Uuid {
public:
    Uuid();
    
    ~Uuid() {}
    
public:
    inline
    String const &
    s_uuid() const { return s_uuid_; }
    
    inline
    bool 
    operator ==(Uuid const & uuid) const {
        return s_uuid_ == uuid.s_uuid_;
    }
 
    inline
    bool
    operator ==(Uuid const & uuid)  {
        return s_uuid_ == uuid.s_uuid_;
    }
    
    inline
    bool
    operator != (Uuid const & uuid) const {
        return s_uuid_ != uuid.s_uuid_;
    }

    
    inline
    bool
    operator != (Uuid & uuid)  {
        return s_uuid_ != uuid.s_uuid_;
    }

    inline
    bool
    operator < (Uuid const & uuid) const {
        return s_uuid_.compare(uuid.s_uuid());
    }
    
private:
    String s_uuid_;
    
};

std::ostream &
operator <<( std::ostream &, Uuid const &);

struct UuidCompare {
    bool operator() (
        Uuid const & u1,
        Uuid const & u2) const {
        return u1.s_uuid() < u2.s_uuid();
    }
};

//std::shared_ptr<Uuid> UuidOP;

#endif
