//
// Created by Joseph Yesselman on 12/31/17.
//

#include <base/json.h>

namespace json {

    JSON JSON::Load( const string &str ) {
        size_t offset = 0;
        return std::move( parse_next( str, offset ) );
    }

}