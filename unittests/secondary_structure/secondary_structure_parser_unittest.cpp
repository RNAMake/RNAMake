

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "secondary_structure/secondary_structure_parser.h"

TEST_CASE( "Test Parsing secondary structure into objects" ) {
    
    SUBCASE("test basic parsing errors") {
        auto p = secondary_structure::Parser();
        
        REQUIRE_THROWS_AS(p.parse("GG+CC", "(((+))"), secondary_structure::Exception);
        REQUIRE_THROWS_AS(p.parse("GG+CC", "()+))"), secondary_structure::Exception);
        REQUIRE_THROWS_AS(p.parse("GGG+CC", "(((+))"), secondary_structure::Exception);
        REQUIRE_THROWS_AS(p.parse("GG+CC", "(.+))"), secondary_structure::Exception);

    }
    
    SUBCASE("test parsing into graph") {
        auto p = secondary_structure::Parser();
        auto g = p.parse("GGG+CC", ".((+))");
        CHECK(g->size() == 5);
        
        p.reset();
        
        SUBCASE("make sure that we are parsing the trailing unpaired residues") {
            g = p.parse("GG+CCAA", "((+))..");
            CHECK(g->size() == 5);
        }
        
    }
    
    SUBCASE("test parsing into single motifs") {
        auto p = secondary_structure::Parser();
        
        auto m = p.parse_to_motif("GG+CC", "((+))");
        CHECK(m->residues().size() == 4);
        CHECK(m->ends().size() == 2);
        CHECK(m->basepairs().size() == 2);
        CHECK(m->chains().size() == 2);
        
        p.reset();
        m = p.parse_to_motif("AGG+CC", ".((+))");
        CHECK(m->ends().size() == 1);
        CHECK(m->basepairs().size() == 2);

        m = p.parse_to_motif("GGAAG+CAAG+CAACC", "((..(+)..(+)..))");
        CHECK(m->ends().size() == 3);
        CHECK(m->basepairs().size() == 4);
        
    }
    
    SUBCASE("test parsing into set of motifs") {
        auto p = secondary_structure::Parser();
        auto motifs = p.parse_to_motifs("GGAGG+CAACCC", "((.((+)..)))");
        CHECK(motifs.size() == 3);
        CHECK(motifs[0]->mtype() == util::MotifType::HELIX);
        CHECK(motifs[1]->mtype() == util::MotifType::TWOWAY);
        
        p.reset();
        motifs = p.parse_to_motifs("GGAAGACAAGACAACC", "((..(.)..(.)..))");
        CHECK(motifs.size() == 4);
        CHECK(motifs[1]->mtype() == util::MotifType::NWAY);

        p.reset();
        motifs = p.parse_to_motifs("GGAAGACAAGACAAGACAAGACCC",
                                   "((..(.)..(.)..(.)..(.)))");
        CHECK(motifs.size() == 6);
        CHECK(motifs[1]->mtype() == util::MotifType::NWAY);


    }
    
    
}