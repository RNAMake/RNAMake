

#include "../common.hpp"

#include <map>

#include <util/uuid.h>

TEST_CASE( "Test unique indentifiers for finding objects" ) {
  using namespace util;
  SUBCASE("test trival") {
    Uuid uuid1 = generate_uuid();
    Uuid uuid2 = generate_uuid();
    CHECK(uuid1 != uuid2);
  }
  SUBCASE("test comparsions") {
    String uuid_str = "3f49aaa8-d40b-003a-2efb-46f2ed9c470b";
    Uuid uuid1 = uuid_from_str(uuid_str);
    Uuid uuid2 = generate_uuid();
    CHECK(uuid1 != uuid2);
    Uuid uuid_copy = uuid1;
    CHECK(uuid_copy.str() == uuid_str);
  }
    
}