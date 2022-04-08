//
// Created by Joe Yesselman on 4/8/22.
//

#include "../common.hpp"

#include <cstdlib>
#include <base/exception.hpp>
#include <base/paths.hpp>

TEST_CASE("test path functions") {
  SUBCASE("test env path") {
    SUBCASE("test trival") {
      CHECK_NOTHROW(base::path::rnamake_path());
    }
    SUBCASE("test not set") {
      String path = base::path::rnamake_path();
      unsetenv("RNAMAKE");
      CHECK_THROWS_AS(base::path::rnamake_path(), base::ResourceException);
      setenv("RNAMAKE", path.c_str(), 1);
    }
    SUBCASE("not a valid path") {
      String path = base::path::rnamake_path();
      setenv("RNAMAKE", "FAKE_PATH", 1);
      CHECK_THROWS_AS(base::path::rnamake_path(), base::ResourceException);
      setenv("RNAMAKE", path.c_str(), 1);
    }
  }
  
}