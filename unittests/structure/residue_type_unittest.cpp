

#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/residue_type_set.h"
#include "structure/is_equal.h"


TEST_CASE("Test ResidueType for Structure" ) {
    
    SUBCASE("Can generate new residue type") {
        auto name = String("GUA");
        auto atom_map = StringIntMap();
        atom_map["P"] = 0;
        auto rt = structure::ResidueType(name, atom_map, structure::SetType::RNA);
        
        CHECK(rt.short_name() == "G");
        CHECK(rt.atom_pos_by_name("P") == 0);
        
    }
    
    SUBCASE("Check matching alternate names") {
        auto name = String("GUA");
        auto atom_map = StringIntMap();
        atom_map["P"] = 0;
        auto rt = structure::ResidueType(name, atom_map, structure::SetType::RNA);
        auto names = Strings{"GUA", "G", "rG"};
        
        for(auto const & n : names) {
            CHECK(rt.match_name(n));
        }
        
        CHECK(rt.match_name("rC") == 0);

    }
    
    auto rts = structure::ResidueTypeSet();
    
    SUBCASE("Does a residue type exist in set") {
        
        SUBCASE("Normal residues should be found by 3-letter name") {
            CHECK(rts.contains_rtype("GUA"));
        }
        
        SUBCASE("short hand names and alt names also should be found") {
            CHECK(rts.contains_rtype("A"));
            CHECK(rts.contains_rtype("rC"));
        }
        
        CHECK(rts.contains_rtype("FAKE") == 0);
        CHECK(rts.contains_rtype("rCC") == 0);
        CHECK(rts.contains_rtype("AA") == 0);

    }
    
    SUBCASE("Getting correct residue type by name from set") {
        auto rt = rts.get_rtype_by_resname("GUA");
        
        CHECK(rt.short_name() == "G");
        REQUIRE_THROWS_AS(rts.get_rtype_by_resname("FAKE"), structure::ResidueTypeException);
        
    }
    
    SUBCASE("can load amino acid residue types") {
        auto rt = rts.get_rtype_by_resname("ARG");
        CHECK(rt.short_name() == "A");
        CHECK(rt.atom_pos_by_name("N") == 0);
    }
    
}
