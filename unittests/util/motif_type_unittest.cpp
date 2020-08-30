
#include "../common.hpp"

#include <base/types.h>
#include <util/motif_type.h>

TEST_CASE( "Test conversions from and to MotifType Enum", "[util::MotifType]" ) {

    SECTION("MotifType => String => MotifType again") {
        for(const auto& motif : {util::MotifType::TWOWAY,
                                util::MotifType::NWAY,
                                util::MotifType::HAIRPIN,
                                util::MotifType::T_T,
                                util::MotifType::T_T_T,
                                util::MotifType::TWOWAY_SEGMENTS,
                                util::MotifType::HELIX,
                                util::MotifType::TCONTACT,
                                util::MotifType::UNKNOWN}) {
            REQUIRE(util::str_to_type(util::type_to_str(motif)) == motif);
        }
    }

    SECTION("String => MotifType => String") {
        for(const auto& str :  {"TWOWAY",
                    "NWAY",
                    "HAIRPIN",
                    "2X_TWOWAY",
                    "3X_TWOWAY",
                    "TWOWAY_SEGMENTS",
                    "HELIX",
                    "TCONTACT",
                    "UNKNOWN"}   ) {
            REQUIRE(str == util::type_to_str(util::str_to_type(str)));
        }
    }

    SECTION("Should throw") {
        for(const auto& bad_str : {"", " NWAY", "not even close"}) {
            REQUIRE_THROWS(
                    util::str_to_type(bad_str)
                    );
        }
    }

}