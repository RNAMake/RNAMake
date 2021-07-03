

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "structure/residue_type_set.h"
#include "util/bead.h"
#include "base/file_io.h"
#include "base/settings.h"
#include "structure/residue.h"
#include "structure/is_equal.h"
#include "util/uuid.h"


TEST_CASE( "Test Residues for Structure" ) {
    auto rts  = structure::ResidueTypeSet();
    auto path = base::unittest_resource_dir() + "residue/test_str_to_residue.dat";
    auto lines= base::get_lines_from_file(path);
    auto residues  = structure::ResidueOPs();
    for(auto const & l : lines) {
        if(l.size() < 10) { break; } // end of file
        auto r = std::make_shared<structure::Residue>(l, rts);
        residues.push_back(r);
    }
    
    SUBCASE("test getting atoms by name") {
        auto r = residues[0];
        auto name = String("C1'");
        auto a = r->get_atom(name);

        CHECK(a != nullptr);
        CHECK(a.get_name() == "C1'");
        
        REQUIRE_THROWS_AS(r->get_atom("fake"), structure::ResidueException);
    }
    
    SUBCASE("are residues detecting connections properly") {
        auto r1 = residues[0];
        auto r2 = residues[1];
        auto r3 = residues[2];
        
        SUBCASE("connecting from 5' to 3'") { CHECK(r1->connected_to(*r2) == 1);  }
        SUBCASE("should not be connected")  { CHECK(r1->connected_to(*r3) == 0);  }
        SUBCASE("connecting from 3' to 5'") { CHECK(r3->connected_to(*r2) == -1); }
        
    }
    
    SUBCASE("are residues generating steric beads properly") {
        auto r = residues[1];
        auto beads = r->beads_;
        
        SUBCASE("produced the right number of beads") { CHECK(beads.size() == 3); }
        
        CHECK(beads[0].get_type() == util::BeadType::PHOS);
        CHECK(beads[1].get_type() == util::BeadType::SUGAR);
        CHECK(beads[2].get_type() == util::BeadType::BASE);

    }
    
    SUBCASE("are residues copying correctly") {
        auto r = residues[0];
        auto r2 = std::make_shared<structure::Residue>(*r);
        
        CHECK(are_residues_equal(r, r2));
    }
    
    SUBCASE("are residues stringifing correctly") {
        auto r = residues[0];
        auto s = r->get_str();
        auto r2 = std::make_shared<structure::Residue>(s, rts);
        
        SUBCASE("residues should be the same but not have the same id") {
            CHECK(structure::are_residues_equal(r, r2, 0));
        }
    }
    
}


