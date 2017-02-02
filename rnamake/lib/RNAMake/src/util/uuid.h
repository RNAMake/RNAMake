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
    bool 
    operator ==(Uuid const & uuid) const {
        return id_ == uuid.id_;
    }
 
    inline
    bool
    operator ==(Uuid const & uuid)  {
        return id_ == uuid.id_;
    }
    
    inline
    bool
    operator != (Uuid const & uuid) const {
        return id_ != uuid.id_;
    }

    
    inline
    bool
    operator != (Uuid & uuid)  {
        return id_ != uuid.id_;
    }

    
private:
    uint64_t id_;
    
};

/*std::ostream &
operator <<( std::ostream &, Uuid const &);

struct UuidCompare {
    bool operator() (
        Uuid const & u1,
        Uuid const & u2) const {
        return u1.s_uuid() < u2.s_uuid();
    }
};*/

typedef std::shared_ptr<Uuid> UuidOP;

#endif
