

#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/residue_type.h"
#include "structure/residue_type_set.h"
#include "structure/is_equal.h"


TEST_CASE("Test ResidueType for Structure" ) {
    
    SUBCASE("Can generate new residue type") {
        auto name = String("GUA");
        auto atom_map = StringIntMap();
        atom_map["P"] = 0;
        auto strings = Strings();
        auto rt = ResidueType(name, atom_map, SetType::RNA, strings);
        
        CHECK(rt.get_short_name() == 'G');
        CHECK(rt.get_atom_index("P") == 0);
        
    }

    //TODO ask joe
    
//    SUBCASE("Check matching alternate names") {
//        auto name = String("GUA");
//        auto atom_map = StringIntMap();
//        atom_map["P"] = 0;
//        auto rt = ResidueType(name, atom_map, SetType::RNA, strings);
//        auto names = Strings{"GUA", "G", "rG"};
//
//        for(auto const & n : names) {
//            CHECK(rt.match_name(n));
//        }
//
//        CHECK(rt.match_name("rC") == 0);
//
//    }
    
    auto rts = structure::ResidueTypeSet();
    
    SUBCASE("Does a residue type exist in set") {
        
        SUBCASE("Normal residues should be found by 3-letter name") {
            CHECK(rts.contains_residue_type("GUA"));
        }
        
        SUBCASE("short hand names and alt names also should be found") {
            CHECK(rts.contains_residue_type("A"));
            CHECK(rts.contains_residue_type("rC"));
        }
        
        CHECK(rts.contains_residue_type("FAKE") == 0);
        CHECK(rts.contains_residue_type("rCC") == 0);
        CHECK(rts.contains_residue_type("AA") == 0);

    }
    
    SUBCASE("Getting correct residue type by name from set") {
        auto rt = rts.get_residue_type("GUA");

        CHECK(rt->get_short_name() == 'G');
        REQUIRE_THROWS_AS(rts.get_residue_type("FAKE"), ResidueTypeException);

    }

    SUBCASE("can load amino acid residue types") {
        auto rt = rts.get_residue_type("ARG");
        CHECK(rt->get_short_name() == 'A');
        CHECK(rt->get_atom_index("N") == 0);
    }
    
}
