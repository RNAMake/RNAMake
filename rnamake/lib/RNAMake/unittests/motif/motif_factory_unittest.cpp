
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
        /*auto path = motif_dirs() + "helices/HELIX.IDEAL";
        auto m = mf.motif_from_file(path);
        
        auto aligned_m = mf.align_motif_to_common_frame(m, 1);
        mf.standardize_motif(aligned_m);
        
        auto dist = aligned_m->ends()[0]->d().magnitude();
        auto r_dist = aligned_m->ends()[0]->r().difference(ref_m->ends()[0]->r());
        
        REQUIRE(dist < 1.0);
        REQUIRE(r_dist < 0.001);
        REQUIRE(m->ends()[1]->uuid() == aligned_m->ends()[0]->uuid());
        REQUIRE(m->sequence() == aligned_m->sequence());
        
        aligned_m = mf.align_motif_to_common_frame(m, 0);
        mf.standardize_motif(aligned_m);
        
        REQUIRE(m->sequence() != aligned_m->sequence());*/
    }
    
    SECTION("test generating motifs from basepairs") {
        auto path = base_dir() + "/rnamake/unittests/resources/motifs/HELIX.IDEAL";
        auto m = mf.motif_from_file(path);
        auto bps = m->basepairs();
        
        auto m_bps = mf.motif_from_bps(bps);
        
        REQUIRE(m_bps->ends().size() == 2);
        REQUIRE(m_bps->sequence() == "GG&CC");
        REQUIRE(m_bps->dot_bracket() == "((&))");
        REQUIRE(m_bps->chains().size() == 2);
        REQUIRE(m_bps->end_ids().size() == 2);
        
        SECTION("will not have two basepair ends should throw error as user is not expecting this") {
        
            bps[1]->bp_type("c...");
            REQUIRE_THROWS_AS(m_bps = mf.motif_from_bps(bps), MotifFactoryException);
        
        }
        
        m_bps = mf.motif_from_bps(BasepairOPs{bps[0]});
        REQUIRE(m_bps->ends().size() == 1);
        
    }
    
    SECTION("Loading in proteins with RNA") {
        auto path = base_dir() + "/rnamake/unittests/resources/pdbs/5g2x.pdb";
        auto m = mf.motif_from_file(path, false, true);
        REQUIRE(m->protein_beads().size() != 0);
    }

    SECTION("test alignment setup") {
        auto path = motif_dirs() + "helices/HELIX.IDEAL.2";
        auto m = mf.motif_from_file(path);

        path = motif_dirs() + "base.motif";
        auto base_motif_1 = file_to_motif(path);
        auto base_motif_2 = file_to_motif(path);

        auto aligned_motif_1 = get_aligned_motif(base_motif_1->ends()[1], m->ends()[0], m);
        auto aligned_motif_2 = get_aligned_motif(aligned_motif_1->ends()[1], base_motif_2->ends()[0], base_motif_2);

        base_motif_1->to_pdb("base.pdb", 1, 1);
        aligned_motif_1->to_pdb("aligned_1.pdb", 1, 1);
        aligned_motif_2->to_pdb("aligned_2.pdb", 1, 1);

    }
    
}
















