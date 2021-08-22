//
// Created by Joseph Yesselman on 1/29/17.
//


#include <primitives/basepair.h>

namespace primitives {

    util::Uuid const &
    Basepair::get_partner(util::Uuid const & r_uuid) const {
        if (r_uuid == _res1_uuid) { return _res2_uuid; }
        else if (r_uuid == _res2_uuid) { return _res1_uuid; }
        else {
            throw BasepairException("called partner with uuid that does not exist in basepair");
        }
    }

}
