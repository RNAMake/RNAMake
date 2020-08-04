


//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "secondary_structure/secondary_structure_parser.h"

TEST_CASE( "Test Parsing secondary structure into objects", "[SSParser]" ) {
    
    SECTION("test basic parsing errors") {
        auto p = secondary_structure::Parser();
        
        REQUIRE_THROWS_AS(p.parse("GG+CC", "(((+))"), secondary_structure::Exception);
        REQUIRE_THROWS_AS(p.parse("GG+CC", "()+))"), secondary_structure::Exception);
        REQUIRE_THROWS_AS(p.parse("GGG+CC", "(((+))"), secondary_structure::Exception);
        REQUIRE_THROWS_AS(p.parse("GG+CC", "(.+))"), secondary_structure::Exception);

    }
    
    SECTION("test parsing into graph") {
        auto p = secondary_structure::Parser();
        auto g = p.parse("GGG+CC", ".((+))");
        REQUIRE(g->size() == 5);
        
        p.reset();
        
        SECTION("make sure that we are parsing the trailing unpaired residues") {
            g = p.parse("GG+CCAA", "((+))..");
            REQUIRE(g->size() == 5);
        }
        
    }
    
    SECTION("test parsing into single motifs") {
        auto p = secondary_structure::Parser();
        
        auto m = p.parse_to_motif("GG+CC", "((+))");
        REQUIRE(m->residues().size() == 4);
        REQUIRE(m->ends().size() == 2);
        REQUIRE(m->basepairs().size() == 2);
        REQUIRE(m->chains().size() == 2);
        
        p.reset();
        m = p.parse_to_motif("AGG+CC", ".((+))");
        REQUIRE(m->ends().size() == 1);
        REQUIRE(m->basepairs().size() == 2);

        m = p.parse_to_motif("GGAAG+CAAG+CAACC", "((..(+)..(+)..))");
        REQUIRE(m->ends().size() == 3);
        REQUIRE(m->basepairs().size() == 4);
        
    }
    
    SECTION("test parsing into set of motifs") {
        auto p = secondary_structure::Parser();
        auto motifs = p.parse_to_motifs("GGAGG+CAACCC", "((.((+)..)))");
        REQUIRE(motifs.size() == 3);
        REQUIRE(motifs[0]->mtype() == util::MotifType::HELIX);
        REQUIRE(motifs[1]->mtype() == util::MotifType::TWOWAY);
        
        p.reset();
        motifs = p.parse_to_motifs("GGAAGACAAGACAACC", "((..(.)..(.)..))");
        REQUIRE(motifs.size() == 4);
        REQUIRE(motifs[1]->mtype() == util::MotifType::NWAY);

        p.reset();
        motifs = p.parse_to_motifs("GGAAGACAAGACAAGACAAGACCC",
                                   "((..(.)..(.)..(.)..(.)))");
        REQUIRE(motifs.size() == 6);
        REQUIRE(motifs[1]->mtype() == util::MotifType::NWAY);


    }
    
    
}