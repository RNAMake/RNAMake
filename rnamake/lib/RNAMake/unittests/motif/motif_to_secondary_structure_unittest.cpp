
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif_to_secondary_structure.h"

TEST_CASE( "Test converting 3D motifs into secondary structure objects", "[MotiftoSS]" ) {
    
    auto path = motif_dirs() + "base.motif";
    auto m = file_to_motif(path);
    path = motif_dirs() + "ref.motif";
    auto ref_m = file_to_motif(path);;
  
    auto converter = MotiftoSecondaryStructure();
    auto ss = converter.to_secondary_structure(m);
    
    REQUIRE(ss->basepairs().size() == m->basepairs().size());
    REQUIRE(ss->residues().size() == m->residues().size());
    
    SECTION("should be able to find all the residues in their secondary structure form ") {
    
        for(auto const & r : m->residues()) {
            REQUIRE(ss->get_residue(r->uuid()) != nullptr);
        }
    
    }
    
    SECTION("should be all to find all the basepairs in their secondary structure form") {
        for(auto const & bp : m->basepairs()) {
            REQUIRE(ss->get_basepair(bp->uuid()).size() == 1);
        }
    }

    ss = converter.to_secondary_structure(ref_m);
    REQUIRE(ss->residues().size() == 2);
    

    
}