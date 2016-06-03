
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/structure.h"
#include "structure/basepair_state.h"
#include "structure/is_equal.hpp"


TEST_CASE( "Test Basepair State for Structure", "[BasepairState]" ) {

    SECTION("test loading basepair from string") {
        auto path = unittest_resource_dir() + "/structure/test_str_to_basepairstate.dat";
        auto lines = get_lines_from_file(path);
        auto bp_state = BasepairStateOP(nullptr);
        
        for(auto const & l : lines) {
            if(l.length() < 5) { break; }
            bp_state = std::make_shared<BasepairState>(l);
        }
     
        auto s = bp_state->to_str();
        auto bp_state_2 = std::make_shared<BasepairState>(s);
    
        REQUIRE(are_basepair_states_equal(*bp_state, *bp_state_2));
        
    }
    
    SECTION("test generating rotation and translation between basepairs") {
        auto path = unittest_resource_dir() + "/structure/get_transforming_r_and_t_test.dat";
        auto lines = get_lines_from_file(path);

        auto dummy = BasepairState();
        int fail = 0;
        for(auto const & l : lines) {
            if(l.length() < 5) { break;}
            auto spl = split_str_by_delimiter(l, "|");
            auto bp_state_1 = BasepairState(spl[0]);
            auto bp_state_2 = BasepairState(spl[1]);
            auto t = Point(spl[2]);
            
            bp_state_1.get_transforming_r_and_t(bp_state_2, dummy);
            auto result_t = t + bp_state_1.d();
            if(! are_xyzVector_equal(result_t, dummy.d())) {
                fail = 1;
            }
            
        }
        
        REQUIRE(fail == 0);
        
    }
    
    
    
}