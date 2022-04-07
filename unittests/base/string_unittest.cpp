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
  SUBCASE("Test splitting") {
    String test = "test1 test2 test3";
    Strings spl = base::string::split(test, " ");
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
    SUBCASE("leading spaces") {
      // trailing space(s)
      const auto raw_line = String{"1 2 3 4 5 6 "};
      const auto target = Strings{"1", "2", "3", "4", "5", "6"};
      const auto actual = base::tokenize_line(raw_line);

      CHECK(target.size() == actual.size());
      auto targ_it = target.cbegin();
      auto actual_it = actual.cbegin();

      for (; targ_it != target.cend(); ++targ_it, ++actual_it) {
        CHECK(*targ_it == *actual_it);
      }
    }

    SUBCASE("spaces on both sides") {
      // leading space(s)
      const auto raw_line = String{"     f s _ $ % "};
      const auto target = Strings{"f", "s", "_", "$", "%"};
      const auto actual = base::tokenize_line(raw_line);

      CHECK(target.size() == actual.size());
      auto targ_it = target.cbegin();
      auto actual_it = actual.cbegin();

      for (; targ_it != target.cend(); ++targ_it, ++actual_it) {
        CHECK(*targ_it == *actual_it);
      }
    }

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
