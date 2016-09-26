
#include "../common.hpp"

#include "data_structure/tree/tree.h"


TEST_CASE( "Test general tree data structure", "[Tree]" ) {

    SECTION("test adding new data for dynamic trees") {
        
        auto t = TreeDynamic<int>();
        t.add_data(0);        //index 0
        t.add_data(1);        //index 1
        t.add_data(2, 0);     //index 2
        
        REQUIRE(t.get_node(0)->children().size() == 2);
        
        SECTION("cannot add to no existent parent") {
            REQUIRE_THROWS_AS(t.add_data(3, 5), TreeException);
        }
        
        t.remove_node(2);
        REQUIRE(t.get_node(0)->children().size() == 1);
        
        SECTION("cannot add to no existent parent") {
            REQUIRE_THROWS_AS(t.add_data(3, 2), TreeException);
        }
        
        t.add_data(2, 0);    //index 3
        
    }
    
    SECTION("check for catching non existant nodes") {
        auto t = TreeDynamic<int>();
        t.add_data(0);        //index 0
        
        auto n = t.get_node(0);
        REQUIRE(n != nullptr);
        REQUIRE_THROWS_AS(t.get_node(2), TreeException);
        REQUIRE_THROWS_AS(t.remove_node(2), TreeException);
        
    }
    
    SECTION("test adding new data for static trees") {
        auto t = TreeStatic<int>();
        t.add_data(0, 2);          //index 0, 2 children
        t.add_data(1, 2, 0, 0);
        t.add_data(2, 2, 0, 1);
        
        SECTION("cannot add node connected to index 0 since it has no free children slots") {
            REQUIRE_THROWS_AS(t.add_data(3, 2, 0), TreeException);
        }
        
        SECTION("cannot add node conected to index 2 does not have a slot 2") {
            REQUIRE_THROWS_AS(t.add_data(3, 2, 2, 2), TreeException);
        }
        
        REQUIRE(t.size() == 3);
        REQUIRE(t.get_node(0)->available_children_pos().size() == 0);
        REQUIRE(t.get_node(1)->available_children_pos().size() == 2);
        REQUIRE(t.get_node(1)->available_pos(1) == 1);

        t.remove_node(2);
        REQUIRE(t.get_node(0)->available_children_pos().size() == 1);
        REQUIRE(t.get_node(0)->available_pos(1) == 1);
    }
    
    SECTION("test remove node level") {
        auto t = TreeStatic<int>();
        t.add_data(0, 2);          //index 0, 2 children
        t.add_data(1, 2, 0, 0);
        t.add_data(2, 2, 0, 1);
        t.increase_level();
        
        t.add_data(3, 2);
        t.add_data(4, 2);
        t.remove_level(1);
        
        REQUIRE(t.size() == 3);
        
        t.add_data(3, 2);
        REQUIRE(t.size() == 4);
        
        t.remove_level(0);
        
        REQUIRE(t.size() == 0);

    
    }
    
    SECTION("test removeing nodes") {
        auto t = TreeStatic<int>();
        t.add_data(0, 2);          //index 0, 2 children
        t.remove_node(0);
        
        REQUIRE(t.last_node() == nullptr);
        
        
    }
 
}