


//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/file_io.h"
#include "base/settings.h"
#include "structure/residue.h"
#include "structure/residue_type_set.h"
#include "structure/is_equal.hpp"


TEST_CASE( "Test Residues for Structure", "[Residue]" ) {
    auto rts  = structure::ResidueTypeSet();
    auto path = base::unittest_resource_dir() + "residue/test_str_to_residue.dat";
    auto lines= base::get_lines_from_file(path);
    auto residues  = structure::ResidueOPs();
    for(auto const & l : lines) {
        if(l.size() < 10) { break; } // end of file
        auto r = std::make_shared<structure::Residue>(l, rts);
        residues.push_back(r);
    }
    
    SECTION("test getting atoms by name") {
        auto r = residues[0];
        auto name = String("C1'");
        auto a = r->get_atom(name);

        REQUIRE(a != nullptr);
        REQUIRE(a->name() == "C1'");
        
        REQUIRE_THROWS_AS(r->get_atom("fake"), structure::ResidueException);
    }
    
    SECTION("are residues detecting connections properly") {
        auto r1 = residues[0];
        auto r2 = residues[1];
        auto r3 = residues[2];
        
        SECTION("connecting from 5' to 3'") { REQUIRE(r1->connected_to(*r2) == 1);  }
        SECTION("should not be connected")  { REQUIRE(r1->connected_to(*r3) == 0);  }
        SECTION("connecting from 3' to 5'") { REQUIRE(r3->connected_to(*r2) == -1); }
        
    }
    
    SECTION("are residues generating steric beads properly") {
        auto r = residues[1];
        auto beads = r->get_beads();
        
        SECTION("produced the right number of beads") { REQUIRE(beads.size() == 3); }
        
        REQUIRE(beads[0].btype() == structure::BeadType::PHOS);
        REQUIRE(beads[1].btype() == structure::BeadType::SUGAR);
        REQUIRE(beads[2].btype() == structure::BeadType::BASE);

    }
    
    SECTION("are residues copying correctly") {
        auto r = residues[0];
        auto r2 = std::make_shared<structure::Residue>(*r);
        
        REQUIRE(are_residues_equal(r, r2));
    }
    
    SECTION("are residues stringifing correctly") {
        auto r = residues[0];
        auto s = r->to_str();
        auto r2 = std::make_shared<structure::Residue>(s, rts);
        
        SECTION("residues should be the same but not have the same id") {
            REQUIRE(structure::are_residues_equal(r, r2, 0));
        }
    }
    
}


