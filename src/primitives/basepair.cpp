//
// Created by Joseph Yesselman on 1/29/17.
//


#include <primitives/basepair.h>

namespace primitives {

    util::Uuid const &
    Basepair::get_partner(util::Uuid const & r_uuid) const {
        if (r_uuid == res1_uuid_) { return res2_uuid_; }
        else if (r_uuid == res2_uuid_) { return res1_uuid_; }
        else {
            throw BasepairException("called partner with uuid that does not exist in basepair");
        }
    }

}
