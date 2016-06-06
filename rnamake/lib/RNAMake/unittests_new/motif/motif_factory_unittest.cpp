
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "motif/motif_factory.h"

TEST_CASE( "Test Motif creation with Motif Factory", "[MotifFactory]" ) {
    
    auto mf = MotifFactory();
    
    auto path = motif_dirs() + "ref.motif";
    auto ref_m = file_to_motif(path);;
    
    SECTION("test loading motif from pdb file") {
        auto path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
        auto m = mf.motif_from_file(path);
        
        REQUIRE(m->residues().size() == 157);
    }
    
    SECTION("load motif from directory") {
        auto path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6";
        auto m = mf.motif_from_file(path);
        
        REQUIRE(m->residues().size() == 157);

    }

    SECTION("returns errors if file does not exist") {
        REQUIRE_THROWS_AS(mf.motif_from_file("fake.pdb"), MotifFactoryException);
        REQUIRE_THROWS_AS(mf.motif_from_file("fake"), MotifFactoryException);
        
        auto path = unittest_resource_dir() + "/motif/empty";
        REQUIRE_THROWS_AS(mf.motif_from_file(path), MotifFactoryException);
    }
 
    SECTION("test standardizing motifs, i.e. making sure they all behave the same") {
        auto path = motif_dirs() + "helices/HELIX.IDEAL";
        auto m = mf.motif_from_file(path);
        
        auto aligned_m = mf.align_motif_to_common_frame(m, 1);
        mf.standardize_motif(aligned_m);
        
        auto dist = aligned_m->ends()[0]->d().magnitude();
        auto r_dist = aligned_m->ends()[0]->r().difference(ref_m->ends()[0]->r());
        
        REQUIRE(dist < 1.0);
        REQUIRE(r_dist < 0.001);
        REQUIRE(m->ends()[1]->uuid() == aligned_m->ends()[0]->uuid());
        //REQUIRE(m->sequence() != aligned_m->sequence());
        
        std::cout << aligned_m->ends()[0]->res1()->name() << std::endl;
        
    }
    
}