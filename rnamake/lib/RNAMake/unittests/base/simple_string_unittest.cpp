
#include "../common.hpp"

#include "base/simple_string.h"

TEST_CASE( "Test very fast and memory efficient string", "[SimpleString]" ) {
    auto test_str = String("test");
    auto ss = SimpleString(test_str);
    auto simple_strings = std::vector<SimpleString>();
    for(int i = 0; i < 10; i++) {
        test_str += std::to_string(i);
        auto ss = SimpleString(test_str);
        simple_strings.push_back(ss);
    }

    int pos = 0;
    for(int i = 0; i < 1000000000; i++) {
        pos = rand() % 10;
        auto ss_2 = SimpleString(simple_strings[pos]);
    }

    //auto ss2 = SimpleString(ss);

}