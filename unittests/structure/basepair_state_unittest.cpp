
#include "../common.hpp"

#include <base/file_io.h>
#include <base/settings.h>
#include <structure/structure.h>
#include <structure/basepair_state.h>
#include <structure/is_equal.h>

TEST_CASE( "Test structure::Basepair State for Structure", "[structure::BasepairState]" ) {

    SECTION("test loading basepair from string") {
        auto path = base::unittest_resource_dir() + "/structure/test_str_to_basepairstate.dat";
        auto lines =base::get_lines_from_file(path);
        auto bp_state = structure::BasepairStateOP(nullptr);
        
        for(auto const & l : lines) {
            if(l.length() < 5) { break; }
            bp_state = std::make_shared<structure::BasepairState>(l);
        }
     
        auto s = bp_state->to_str();
        auto bp_state_2 = std::make_shared<structure::BasepairState>(s);
    
        REQUIRE(structure::are_basepair_states_equal(*bp_state, *bp_state_2));
        
    }
    
    SECTION("test generating rotation and translation between basepairs") {
        auto path = base::unittest_resource_dir() + "/structure/get_transforming_r_and_t_test.dat";
        auto lines =base::get_lines_from_file(path);

        auto dummy = structure::BasepairState();
        int fail = 0;
        int i = -1;
        for(auto const & l : lines) {
            i++;
            if(l.length() < 5) { break;}
            auto spl = base::split_str_by_delimiter(l, "|");
            auto bp_state_1 = structure::BasepairState(spl[0]);
            auto bp_state_2 = structure::BasepairState(spl[1]);
            auto t = math::Point(spl[2]);
            auto r = math::Matrix(spl[3]);
            
            bp_state_1.get_transforming_r_and_t(bp_state_2, dummy);
            auto result_t = t + bp_state_1.d();
            
            if(! math::are_xyzVector_equal(result_t, dummy.d())) {
                fail = 1;
            }
            
            if(! math::are_xyzMatrix_equal(r, dummy.r())) {
                fail = 1;
            }
            
        }
        
        REQUIRE(fail == 0);
        
    }
    
    SECTION("test random tranformation between basepairs") {
        auto path = base::unittest_resource_dir() + "/structure/test_get_transformed_state.dat";
        auto lines =base::get_lines_from_file(path);
        auto dummy = structure::BasepairState();
        auto dummy_2 = structure::BasepairState();
        
        int fail = 0;
        
        for(auto const & l : lines) {
            if(l.length() < 5) { break;}
            auto spl = base::split_str_by_delimiter(l, "|");
            auto bp_state_1 = structure::BasepairState(spl[0]);
            auto bp_state_2 = structure::BasepairState(spl[1]);
            auto bp_result  = structure::BasepairState(spl[2]);
            
            bp_state_1.get_transforming_r_and_t(bp_state_2, dummy);
            bp_state_2.get_transformed_state(dummy, dummy_2);

            if(! structure::are_basepair_states_equal(bp_state_1, dummy_2)) {
                fail = 1;
            }
            if(! structure::are_basepair_states_equal(bp_result, dummy_2)) {
                fail = 1;
            }
        }
        
        REQUIRE(fail == 0);

    }
    
    SECTION("test repeated alignment") {
        auto path = base::unittest_resource_dir() + "/structure/test_bp_state_align_1.dat";
        auto lines =base::get_lines_from_file(path);
        
        auto spl = base::split_str_by_delimiter(lines[0], "|");
        auto bp_state_1 = structure::BasepairState(spl[0]);
        auto bp_state_2 = structure::BasepairState(spl[1]);
        auto dummy = structure::BasepairState();
        auto dummy_2 = structure::BasepairState();
        
        auto dist = 0.0;
        auto rdist = 0.0;
        int fail = 0;
        
        for(int i = 0; i < 100; i++) {
            bp_state_1.get_transforming_r_and_t(bp_state_2, dummy);
            bp_state_2.get_transformed_state(dummy, dummy_2);
           
            bp_state_2.set(dummy_2);
            
            dist = bp_state_1.d().distance(bp_state_2.d());
            rdist = bp_state_1.r().difference(bp_state_2.r());
            
            if(dist > 0.001) { fail =1; }
            if(rdist > 0.001) { fail = 1; }
            
            //std::cout << dist << " " << rdist << std::endl;
            
        }
        
        REQUIRE(fail == 0);
        fail = 0;
        
        bp_state_1 = structure::BasepairState(spl[1]);
        bp_state_2 = structure::BasepairState(spl[0]);

        for(int i = 0; i < 100; i++) {
            bp_state_1.get_transforming_r_and_t(bp_state_2, dummy);
            bp_state_2.get_transformed_state(dummy, dummy_2);
            
            bp_state_2.set(dummy_2);
            
            dist = bp_state_1.d().distance(bp_state_2.d());
            rdist = bp_state_1.r().difference(bp_state_2.r());
            
            if(dist > 0.001) { fail =1; }
            if(rdist > 0.001) { fail = 1; }
            
            //std::cout << dist << " " << rdist << std::endl;
            
        }

        REQUIRE(fail == 0);

        
    }
    
    SECTION("test transform() instead of get_transformed_state()") {
        auto path = base::unittest_resource_dir() + "/structure/test_bp_state_align_1.dat";
        auto lines =base::get_lines_from_file(path);

        auto dist = 0.0;
        auto rdist = 0.0;
        int fail = 0;

        for(int i = 0; i < 100; i++) {
            auto spl = base::split_str_by_delimiter(lines[0], "|");
            auto bp_state_1 = structure::BasepairState(spl[0]);
            auto bp_state_2 = structure::BasepairState(spl[1]);

            auto dummy = structure::BasepairState();
            auto dummy_2 = structure::BasepairState();

            bp_state_1.get_transforming_r_and_t(bp_state_2, dummy);
            bp_state_2.get_transformed_state(dummy, dummy_2);

            auto bp_state_copy = bp_state_2.copy();
            bp_state_2.set(dummy_2);

            bp_state_copy.transform(dummy.r().transposed(), dummy.d());


            dist = bp_state_copy.d().distance(bp_state_2.d());
            rdist = bp_state_copy.r().difference(bp_state_2.r());

            if(dist > 0.001) { fail =1; }
            if(rdist > 0.001) { fail = 1; }

            REQUIRE(!fail);

            //std::cout << dist << " " << rdist << std::endl;

        }
    }


}




































