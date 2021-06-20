
#include "../common.hpp"

#include "data_structure/tree/tree.h"

TEST_CASE( "Test general tree data structure") {

    SUBCASE("test adding new data for dynamic trees") {
        auto t = data_structure::tree::TreeDynamic<int>();
        t.add_data(0);        //index 0
        t.add_data(1);        //index 1
        t.add_data(2, 0);     //index 2
        
        CHECK(t.get_node(0)->children().size() == 2);
        
        SUBCASE("cannot add to no existent parent") {
            REQUIRE_THROWS_AS(t.add_data(3, 5), data_structure::tree::TreeException);
        }
        
        t.remove_node(2);
        CHECK(t.get_node(0)->children().size() == 1);
        
        SUBCASE("cannot add to no existent parent") {
            REQUIRE_THROWS_AS(t.add_data(3, 2), data_structure::tree::TreeException);
        }
        
        t.add_data(2, 0);    //index 3
        
    }
    
    SUBCASE("check for catching non existant nodes") {
        auto t = data_structure::tree::TreeDynamic<int>();
        t.add_data(0);        //index 0
        
        auto n = t.get_node(0);
        CHECK(n != nullptr);
        REQUIRE_THROWS_AS(t.get_node(2), data_structure::tree::TreeException);
        REQUIRE_THROWS_AS(t.remove_node(2), data_structure::tree::TreeException);
        
    }
    
    SUBCASE("test adding new data for static trees") {
        auto t = data_structure::tree::TreeStatic<int>();
        t.add_data(0, 2);          //index 0, 2 children
        t.add_data(1, 2, 0, 0);
        t.add_data(2, 2, 0, 1);
        
        SUBCASE("cannot add node connected to index 0 since it has no free children slots") {
            REQUIRE_THROWS_AS(t.add_data(3, 2, 0), data_structure::tree::TreeException);
        }
        
        SUBCASE("cannot add node conected to index 2 does not have a slot 2") {
            REQUIRE_THROWS_AS(t.add_data(3, 2, 2, 2), data_structure::tree::TreeException);
        }
        
        CHECK(t.size() == 3);
        CHECK(t.get_node(0)->available_children_pos().size() == 0);
        CHECK(t.get_node(1)->available_children_pos().size() == 2);
        CHECK(t.get_node(1)->available_pos(1) == 1);

        t.remove_node(2);
        CHECK(t.get_node(0)->available_children_pos().size() == 1);
        CHECK(t.get_node(0)->available_pos(1) == 1);
    }
    
    SUBCASE("test remove node level") {
        auto t = data_structure::tree::TreeStatic<int>();
        t.add_data(0, 2);          //index 0, 2 children
        t.add_data(1, 2, 0, 0);
        t.add_data(2, 2, 0, 1);
        t.increase_level();
        
        t.add_data(3, 2);
        t.add_data(4, 2);
        t.remove_level(1);
        
        CHECK(t.size() == 3);
        
        t.add_data(3, 2);
        CHECK(t.size() == 4);
        
        t.remove_level(0);
        
        CHECK(t.size() == 0);

    
    }
    
    SUBCASE("test removeing nodes") {
        auto t = data_structure::tree::TreeStatic<int>();
        t.add_data(0, 2);          //index 0, 2 children
        t.remove_node(0);
        
        CHECK(t.last_node() == nullptr);
        
        
    }
 
}