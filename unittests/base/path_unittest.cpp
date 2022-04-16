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
    SUBCASE("test resources path") {
      String path = base::path::resources_path();
      CHECK(path == base::path::rnamake_path() + "resources/");
    }
  }
  SUBCASE("test filename") {
    SUBCASE("test valid filename") {
      String path_str = base::path::unittest_resource_path() + "base/path";
      String filename = base::path::filename(path_str);
      CHECK(filename == "path");
    }
    SUBCASE("empty string") {
      String path_str;
      CHECK_THROWS_AS(base::path::filename(path_str), base::InputException);
    }
    SUBCASE("has no file name") {
      String path_str = base::path::unittest_resource_path() + "base/path/";
      CHECK_THROWS_AS(base::path::filename(path_str), base::InputException);
    }
  }
  SUBCASE("test parent_path") {
    SUBCASE("test valid parent_path") {
      String path_str = base::path::unittest_resource_path() + "base/path";
      String base = base::path::unittest_resource_path() + "base";
      String parent_path = base::path::parent_dir(path_str);
      CHECK(parent_path == base);
    }
    SUBCASE("empty string") {
      String path_str;
      CHECK_THROWS_AS(base::path::parent_dir(path_str), base::InputException);
    }
    SUBCASE("has no parent_path") {
      String path_str = "/";
      CHECK(base::path::parent_dir(path_str) == "/");
    }
  }
  SUBCASE("test get file contents") {
    SUBCASE("test valid file") {
      String path_str = base::path::unittest_resource_path() +
                        "base/path/test.csv";
      Strings lines;
      base::path::get_lines_from_file(path_str, lines);
      std::cout << lines.size() << std::endl;

    }
  }
  
}