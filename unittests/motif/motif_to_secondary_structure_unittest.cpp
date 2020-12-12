

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif_to_secondary_structure.h"

#include <util/find_pair.h> // Is this import necessary?

TEST_CASE( "Test converting 3D motifs into secondary structure objects" ) {
    
    auto path = base::motif_dirs() + "base.motif";
    auto m = motif::file_to_motif(path);
    path = base::motif_dirs() + "ref.motif";
    auto ref_m = motif::file_to_motif(path);;
  
    auto converter = motif::MotiftoSecondaryStructure();
    auto ss = converter.to_secondary_structure(m);
    
    CHECK(ss->basepairs().size() == m->basepairs().size());
    CHECK(ss->residues().size() == m->residues().size());
    
    SUBCASE("should be able to find all the residues in their secondary structure form ") {
    
        for(auto const & r : m->residues()) {
            CHECK(ss->get_residue(r->uuid()) != nullptr);
        }
    
    }
    
    SUBCASE("should be all to find all the basepairs in their secondary structure form") {
        for(auto const & bp : m->basepairs()) {
            CHECK(ss->get_basepair(bp->uuid()).size() == 1);
        }
    }

    ss = converter.to_secondary_structure(ref_m);
    CHECK(ss->residues().size() == 2);
    

    
}