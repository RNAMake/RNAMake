
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/residue_type_set.h"
#include "structure/is_equal.h"


TEST_CASE("Test ResidueType for Structure", "[ResidueType]" ) {
    
    SECTION("Can generate new residue type") {
        auto name = String("GUA");
        auto atom_map = StringIntMap();
        atom_map["P"] = 0;
        auto rt = structure::ResidueType(name, atom_map, structure::SetType::RNA);
        
        REQUIRE(rt.short_name() == "G");
        REQUIRE(rt.atom_pos_by_name("P") == 0);
        
    }
    
    SECTION("Check matching alternate names") {
        auto name = String("GUA");
        auto atom_map = StringIntMap();
        atom_map["P"] = 0;
        auto rt = structure::ResidueType(name, atom_map, structure::SetType::RNA);
        auto names = Strings{"GUA", "G", "rG"};
        
        for(auto const & n : names) {
            REQUIRE(rt.match_name(n));
        }
        
        REQUIRE(rt.match_name("rC") == 0);

    }
    
    auto rts = structure::ResidueTypeSet();
    
    SECTION("Does a residue type exist in set") {
        
        SECTION("Normal residues should be found by 3-letter name") {
            REQUIRE(rts.contains_rtype("GUA"));
        }
        
        SECTION("short hand names and alt names also should be found") {
            REQUIRE(rts.contains_rtype("A"));
            REQUIRE(rts.contains_rtype("rC"));
        }
        
        REQUIRE(rts.contains_rtype("FAKE") == 0);
        REQUIRE(rts.contains_rtype("rCC") == 0);
        REQUIRE(rts.contains_rtype("AA") == 0);

    }
    
    SECTION("Getting correct residue type by name from set") {
        auto rt = rts.get_rtype_by_resname("GUA");
        
        REQUIRE(rt.short_name() == "G");
        REQUIRE_THROWS_AS(rts.get_rtype_by_resname("FAKE"), structure::ResidueTypeException);
        
    }
    
    SECTION("can load amino acid residue types") {
        auto rt = rts.get_rtype_by_resname("ARG");
        REQUIRE(rt.short_name() == "A");
        REQUIRE(rt.atom_pos_by_name("N") == 0);
    }
    
}
