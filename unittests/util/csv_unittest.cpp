//
// Created by Joseph Yesselman on 2/16/20.
//

#include "../common.hpp"

#include <base/paths.hpp>
#include <util/csv.h>

TEST_CASE("Test csv parsing") {
  SUBCASE("test reading a csv") {
    auto path = base::path::unittest_resource_path() + "/base/test.csv";
    auto csv_reader = util::CSVReader<3>(path);
    String test1, test2, test3;
    int count = 0;
    while (csv_reader.read_row(test1, test2, test3)) {
      count += 1;
    }
    CHECK(count == 3);
  }
}