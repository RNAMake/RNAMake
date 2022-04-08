//
// Created by Joseph Yesselman on 2/16/20.
//

//
// Created by Joseph Yesselman on 3/10/19.
//

#include <base/string.hpp>
#include <base/types.hpp>
#include <iostream>

#include "../common.hpp"

TEST_CASE("test string functions") {
  SUBCASE("test string splitting") {
    SUBCASE("Test splitting simple") {
      String test = "test1 test2 test3";
      Strings spl = base::string::split(test, " ");
      CHECK(spl.size() == 3);
      CHECK(spl[0] == "test1");
      CHECK(spl[1] == "test2");
      CHECK(spl[2] == "test3");
    }
    SUBCASE("Two delimiters in a row ") {
      String test = "test1  test2  test3";
      Strings spl = base::string::split(test, " ");
      CHECK(spl.size() == 5);
      spl = base::string::split(test, "  ");
      CHECK(spl.size() == 3);
      CHECK(spl[0] == "test1");
    }
    SUBCASE("Only deliminters") {
      String test = "  ";
      Strings spl = base::string::split(test, " ");
      CHECK(spl.size() == 2);
      CHECK(spl[0] == "");
    }
    SUBCASE("No delimiter") {
      String test = "test1  test2  test3";
      Strings spl = base::string::split(test, "\t");
      CHECK(spl.size() == 1);
      CHECK(spl[0] == test);
    }
    SUBCASE("Empty string") {
      String test;
      Strings spl = base::string::split(test, " ");
      CHECK(spl.size() == 1);
    }
    SUBCASE("escape characters") {
      String test = R"("\0\t\t\t\a")";
      Strings spl = base::string::split(test, R"("\t")");
      CHECK(spl.size() == 1);
    }
  }
  SUBCASE("test string joining") {
    SUBCASE("test trivial case") {
      Strings spl = {"test1", "test2", "test3"};
      String test = base::string::join(spl, " ");
      CHECK(test == "test1 test2 test3");
    }
    SUBCASE("test empty vector") {
      Strings spl = {};
      String test = base::string::join(spl, " ");
      CHECK(test == "");
    }
    SUBCASE("test line joining") {
      Strings spl = {"test1", "test2", "test3"};
      String test = base::string::join(spl, "\n");
      CHECK(test == "test1\ntest2\ntest3");
    }
  }
  SUBCASE("Testing trim methods") {
    SUBCASE("left trim") {
      const String target = {"trimmed"};
      String untrimmed = {" trimmed"};
      CHECK(base::string::left_trim(untrimmed) == target);
      untrimmed = "      trimmed";
      CHECK(base::string::left_trim(untrimmed) == target);
      untrimmed = "\t\ttrimmed";
      CHECK(base::string::left_trim(untrimmed) == target);
      untrimmed = "\n\ttrimmed";
      CHECK(base::string::left_trim(untrimmed) == target);
      String empty_str = {"            "};
      CHECK(base::string::left_trim(empty_str).empty());
    }
    SUBCASE("right trim") {
      const String target = {"trimmed"};
      String untrimmed = {"trimmed "};
      CHECK(base::string::right_trim(untrimmed) == target);
      untrimmed = "trimmed            ";
      CHECK(base::string::right_trim(untrimmed) == target);
      untrimmed = "trimmed\t\t";
      CHECK(base::string::right_trim(untrimmed) == target);
      untrimmed = "trimmed\n\t";
      CHECK(base::string::right_trim(untrimmed) == target);
      String empty_str = {"            "};
      CHECK(base::string::right_trim(empty_str).empty());
    }
    SUBCASE("trim") {
      const String target = {"trimmed"};
      String untrimmed = {" trimmed "};
      CHECK(base::string::trim(untrimmed) == target);
      untrimmed = "         trimmed            ";
      CHECK(base::string::trim(untrimmed) == target);
      untrimmed = "\t\ttrimmed\t\t";
      CHECK(base::string::trim(untrimmed) == target);
      untrimmed = "\n\ttrimmed\n\t";
      CHECK(base::string::trim(untrimmed) == target);
      String empty_str = {"            "};
      CHECK(base::string::trim(empty_str).empty());
    }
  }
  SUBCASE("quote behavior") {
    SUBCASE("test default behavior") {
      const String target = "\'test\'";
      String initial = "test";
      CHECK(base::string::quoted(initial) == target);
    }
    SUBCASE("empty string") {
      const String target ="\'\'";
      String initial;
      CHECK(base::string::quoted(initial) == target);
    }
  }

}
