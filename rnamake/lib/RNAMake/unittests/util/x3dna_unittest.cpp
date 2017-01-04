
#include "../common.hpp"

#include "base/settings.h"
#include "base/file_io.h"
#include "util/x3dna.h"

TEST_CASE( "Test wrapper for x3dna calls", "[X3DNA]" ) {
   
    SECTION("test generating ref_frame file from pdb") {
    
        auto m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6_2/p4p6_2.pdb";
        auto x = X3dna();
        x.generate_ref_frame(m_path);
    
        REQUIRE(file_exists("ref_frames.dat"));
        REQUIRE(!file_exists("basepairs.pdb"));
        
        std::remove("ref_frames.dat");
        
        SECTION("throw error if the pdb does not exist") {
            REQUIRE_THROWS_AS(x.generate_ref_frame("fake.pdb"), X3dnaException);
            
        }
        
        
    }
    
    SECTION("test generating dssr file from pdb") {
        auto m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6_2/p4p6_2.pdb";
        auto x = X3dna();
        x.generate_dssr_file(m_path);
        
        REQUIRE(file_exists("p4p6_2_dssr.out"));
        
        std::remove("p4p6_2_dssr.out");
 
    }
    
    SECTION("test x3dna residue equality") {
        auto res1 = X3Residue(1, "A", "");
        auto res2 = X3Residue(1, "A", "");
        auto res3 = X3Residue(2, "A", "");
        auto res4 = X3Residue(1, "B", "");
        
        REQUIRE(res1 == res2);
        REQUIRE(!(res1 == res3));
        REQUIRE(!(res1 == res4));
    }
    
    SECTION("test getting x3dna basepairs from pdb") {
        
        auto x = X3dna();
        SECTION("should not build new ref frame and dssr files since they already exist") {
        
            auto m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
            auto basepairs = x.get_basepairs(m_path);
        
            REQUIRE(!file_exists("ref_frames.dat"));
            REQUIRE(!file_exists("p4p6_dssr.out"));
        }
        
        SECTION("needs to build new ref frame and dssr files") {
            
            auto m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6_2/p4p6_2.pdb";
            auto basepairs = x.get_basepairs(m_path);
            
            REQUIRE(file_exists("ref_frames.dat"));
            REQUIRE(file_exists("p4p6_2_dssr.out"));
            
            std::remove("ref_frames.dat");
            std::remove("p4p6_2_dssr.out");

        }
        
        SECTION("force creation fo ref_frames file") {
            auto m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
            auto basepairs = x.get_basepairs(m_path, true);
            
            REQUIRE(file_exists("ref_frames.dat"));
            REQUIRE(file_exists("p4p6_dssr.out"));
            
            std::remove("ref_frames.dat");
            std::remove("p4p6_dssr.out");

        }
    }
    
    SECTION("test getting motifs from large rna structure") {
        auto m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
        auto x = X3dna();
        auto motifs = x.get_motifs(m_path);
        REQUIRE(motifs.size() == 22);
    }
    
}
