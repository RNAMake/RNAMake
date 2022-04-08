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
  }

  /*SUBCASE("Testing trim methods") {
    SUBCASE("left trim") {
      const auto target = String{"trimmed"};
      auto untrimmed = String{" trimmed"};
      CHECK(base::ltrim(untrimmed) == target);
      untrimmed = "      trimmed";
      CHECK(base::ltrim(untrimmed) == target);
      untrimmed = "\t\ttrimmed";
      CHECK(base::ltrim(untrimmed) == target);
      untrimmed = "\n\ttrimmed";
      CHECK(base::ltrim(untrimmed) == target);
      auto empty_str = String{"            "};
      CHECK(base::ltrim(empty_str).empty());
    }

    SUBCASE("right trim") {
      const auto target = String{"trimmed"};
      auto untrimmed = String{"trimmed "};
      CHECK(base::rtrim(untrimmed) == target);
      untrimmed = "trimmed            ";
      CHECK(base::rtrim(untrimmed) == target);
      untrimmed = "trimmed\t\t";
      CHECK(base::rtrim(untrimmed) == target);
      untrimmed = "trimmed\n\t";
      CHECK(base::rtrim(untrimmed) == target);
      auto empty_str = String{"            "};
      CHECK(base::rtrim(empty_str).empty());
    }

    SUBCASE("trim") {
      const auto target("trimmed");
      auto untrimmed = String{" trimmed "};
      CHECK(base::trim(untrimmed) == target);
      untrimmed = "         trimmed            ";
      CHECK(base::trim(untrimmed) == target);
      untrimmed = "\t\ttrimmed\t\t";
      CHECK(base::trim(untrimmed) == target);
      untrimmed = "\n\ttrimmed\n\t";
      CHECK(base::trim(untrimmed) == target);
      auto empty_str = String{"            "};
      CHECK(base::trim(empty_str).empty());
    }
  }

  SUBCASE("testing string splitting methods") {

    SUBCASE("empty string") {
      // leading space(s)
      const auto raw_line = String{" "};
      const auto target = Strings{};
      const auto actual = base::tokenize_line(raw_line);

      CHECK(target.size() == actual.size());
      auto targ_it = target.cbegin();
      auto actual_it = actual.cbegin();

      for (; targ_it != target.cend(); ++targ_it, ++actual_it) {
        CHECK(*targ_it == *actual_it);
      }
    }

    SUBCASE("escape characters") {
      const auto raw_line = String{"\0\t\t\t\a"};
      const auto target = Strings{};
      const auto actual = base::tokenize_line(raw_line);

      CHECK(target.size() == actual.size());
      auto targ_it = target.cbegin();
      auto actual_it = actual.cbegin();

      for (; targ_it != target.cend(); ++targ_it, ++actual_it) {
        CHECK(*targ_it == *actual_it);
      }
    }

    SUBCASE("quote behavior") {
      const auto raw_line = String{"\'this has spaces\' last_token"};
      const auto target = Strings{"this has spaces", "last_token"};
      const auto actual = base::tokenize_line(raw_line);

      CHECK(target.size() == actual.size());
      auto targ_it = target.cbegin();
      auto actual_it = actual.cbegin();

      for (; targ_it != target.cend(); ++targ_it, ++actual_it) {
        CHECK(*targ_it == *actual_it);
      }
    }
  }*/
}
