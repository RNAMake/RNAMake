//
// Created by Joseph Yesselman on 2/16/20.
//


#include "../common.hpp"
#include <nlohmann/json.hpp>
#include <limits>


#include <base/settings.h>
#include <base/sys_interface.h>
#include <math/numerical.h>
#include <math/quaternion.h>
#include <util/dssr.h>

TEST_CASE("get_TYPE methods for nlohman::json parsing") {

  SUBCASE("Primitive data types, DOES have key") {
        nlohmann::json json = {
            {"string","this is a string"},
            {"double", 1996.10},
            {"char", "H"},
            {"int", 9000}

        };
        // testing pulls 
        CHECK(util::get_string(json,"string") == "this is a string");
        CHECK(util::get_string(json,"string") != "thisis a string");
            
        CHECK(math::roughly_equal(util::get_double(json,"double"),1996.10));
        CHECK(!math::roughly_equal(util::get_double(json,"double"),2000.));
    
        CHECK(util::get_char(json,"char") == 'H');
        CHECK(util::get_char(json,"char") != 'h');

        CHECK(util::get_int(json,"int") == 9000);
        CHECK(util::get_int(json,"int") != 900);
    }

    SUBCASE("Primitive data types, DOESN'T have key") {
        nlohmann::json json = {};
        
        CHECK(util::get_string(json,"not real") == "NA");
        
        CHECK(math::roughly_equal(util::get_double(json,"not real"),-1.));
        
        CHECK(util::get_char(json,"not real") == ' ');
        
        CHECK(util::get_int(json, "not real") == std::numeric_limits<int>::min());
        
    }

    SUBCASE("Primitive data types, null values") {
        nlohmann::json json = "{\"string\" : null, \"double\" : null, \"int\" : null, \"char\" : null}"_json;
        
        CHECK(util::get_string(json,"string")  == "NA");        
        CHECK(util::get_int(json,"int")  == std::numeric_limits<int>::min());        
        CHECK(util::get_char(json,"char")  == ' ');        
        CHECK(math::roughly_equal(util::get_double(json,"double"),-1.));        

    }
//TODO composite datatypes
//math::Point
//get_point(const nlohmann::json& , const String&);
//
//math::Quaternion
//get_quaternion(const nlohmann::json&, const String&);
//
//math::Matrix
//get_matrix(const nlohmann::json&);
//
//Reals
//get_reals(const nlohmann::json&, const String&);
//
//Ints
//get_ints(const nlohmann::json&, const String&);

}

