
#include "../common.hpp"

#include <map>

#include "base/settings.h"
#include "base/file_io.h"
#include "util/uuid.h"

TEST_CASE( "Test unique indentifiers for finding objects", "[Uuid]" ) {
    
    SECTION("test comparing unique indentifiers") {
        auto u1 = Uuid(), u2 = Uuid();
        auto u3 = u1;
        
        REQUIRE(!(u1 == u2));
        REQUIRE(u1 == u3);
    }
    
    SECTION("test using uuids with maps") {
        auto uuid_int_map = std::map<Uuid, int, UuidCompare>();
        auto u1 = Uuid();
        uuid_int_map[u1] = 1;
        
        for(int i = 0; i < 1000; i++) {
            auto u = Uuid();
            uuid_int_map[u] = 1;
        }
        
        REQUIRE(uuid_int_map.find(u1) != uuid_int_map.end());
    }

}