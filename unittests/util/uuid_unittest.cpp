

#include "../common.hpp"

#include <map>

#include "base/settings.h"
#include "base/file_io.h"
#include "util/uuid.h"

TEST_CASE( "Test unique indentifiers for finding objects" ) {
    
    SUBCASE("test comparing unique indentifiers") {
        auto u1 = util::Uuid(), u2 = util::Uuid();
        auto u3 = u1;
        
        CHECK(!(u1 == u2));
        CHECK(u1 == u3);
    }
    
    SUBCASE("test using uuids with maps") {
        auto uuid_int_map = std::map<util::Uuid, int, util::UuidCompare>();
        auto u1 = util::Uuid();
        uuid_int_map[u1] = 1;
        
        for(int i = 0; i < 1000; i++) {
            auto u = util::Uuid();
            uuid_int_map[u] = 1;
        }
        
        CHECK(uuid_int_map.find(u1) != uuid_int_map.end());
    }

}