

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "motif/motif_factory.h"

TEST_CASE( "Test Motif creation with Motif Factory" ) {
    
    auto mf = motif::MotifFactory();
    
    auto path = base::motif_dirs() + "ref.motif";
    auto ref_m = motif::file_to_motif(path);;
    
    SUBCASE("test loading motif from pdb file") {
        auto path = base::base_dir() + "/unittests/unittest_resources/motifs/p4p6/p4p6.pdb";
        auto m = mf.motif_from_file(path);
        
        CHECK(m->residues().size() == 157);
    }
    
    SUBCASE("load motif from directory") {
        auto path = base::base_dir() + "/unittests/unittest_resources/motifs/p4p6";
        auto m = mf.motif_from_file(path);
        
        CHECK(m->residues().size() == 157);

    }

    SUBCASE("returns errors if file does not exist") {
        REQUIRE_THROWS_AS(mf.motif_from_file("fake.pdb"), motif::MotifFactoryException);
        REQUIRE_THROWS_AS(mf.motif_from_file("fake"), motif::MotifFactoryException);
        
        auto path = base::unittest_resource_dir() + "/motif/empty";
        REQUIRE_THROWS_AS(mf.motif_from_file(path), motif::MotifFactoryException);
    }
 
    SUBCASE("test standardizing motifs, i.e. making sure they all behave the same") {
        /*auto path = base::motif_dirs() + "helices/HELIX.IDEAL";
        auto m = mf.motif_from_file(path);
        
        auto aligned_m = mf.align_motif_to_common_frame(m, 1);
        mf.standardize_motif(aligned_m);
        
        auto dist = aligned_m->ends()[0]->d().magnitude();
        auto r_dist = aligned_m->ends()[0]->r().difference(ref_m->ends()[0]->r());
        
        CHECK(dist < 1.0);
        CHECK(r_dist < 0.001);
        CHECK(m->ends()[1]->uuid() == aligned_m->ends()[0]->uuid());
        CHECK(m->sequence() == aligned_m->sequence());
        
        aligned_m = mf.align_motif_to_common_frame(m, 0);
        mf.standardize_motif(aligned_m);
        
        CHECK(m->sequence() != aligned_m->sequence());*/
    }
    
    SUBCASE("test generating motifs from basepairs") {
        auto path = base::base_dir() + "/unittests/unittest_resources/motifs/HELIX.IDEAL";
        auto m = mf.motif_from_file(path);
        auto bps = m->basepairs();
        
        auto m_bps = mf.motif_from_bps(bps);
        
        CHECK(m_bps->ends().size() == 2);
        CHECK(m_bps->sequence() == "GG&CC");
        CHECK(m_bps->dot_bracket() == "((&))");
        CHECK(m_bps->chains().size() == 2);
        CHECK(m_bps->end_ids().size() == 2);
        
        SUBCASE("will not have two basepair ends should throw error as user is not expecting this") {
        
            bps[1]->bp_type("c...");
            REQUIRE_THROWS_AS(m_bps = mf.motif_from_bps(bps), motif::MotifFactoryException);
        
        }
        
        m_bps = mf.motif_from_bps(structure::BasepairOPs{bps[0]});
        CHECK(m_bps->ends().size() == 1);
        
    }
    
    SUBCASE("Loading in proteins with RNA") {
        auto path = base::base_dir() + "/unittests/unittest_resources/pdbs/5g2x.pdb";
        auto m = mf.motif_from_file(path, false, true);
        CHECK(m->protein_beads().size() != 0);
    }

    SUBCASE("test alignment setup") {
        auto path = base::base_dir() + "/unittests/unittest_resources/motifs/HELIX.IDEAL";//.2";
        auto m = mf.motif_from_file(path);

        path = base::motif_dirs() + "base.motif";
        auto base_motif_1 = motif::file_to_motif(path);
        auto base_motif_2 = motif::file_to_motif(path);

        auto aligned_motif_1 = get_aligned_motif(base_motif_1->ends()[1], m->ends()[0], m);
        auto aligned_motif_2 = get_aligned_motif(aligned_motif_1->ends()[1], base_motif_2->ends()[0], base_motif_2);

        base_motif_1->to_pdb("base.pdb", 1, 1);
        aligned_motif_1->to_pdb("aligned_1.pdb", 1, 1);
        aligned_motif_2->to_pdb("aligned_2.pdb", 1, 1);

    }

    SUBCASE("test forcing set number of chains") {
        auto path = base::base_dir() + "/unittests/unittest_resources/motif/construct_3.pdb";
        auto m = mf.motif_from_file(path, false, true, 1);

        CHECK(m->chains().size() == 1);

    }

}
















