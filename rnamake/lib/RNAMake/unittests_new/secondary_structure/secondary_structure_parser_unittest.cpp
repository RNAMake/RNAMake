


//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "util/uuid.h"
#include "secondary_structure/secondary_structure_parser.h"

TEST_CASE( "Test Parsing secondary structure into objects", "[SSParser]" ) {
    
    SECTION("test basic parsing errors") {
        auto p = sstruct::SecondaryStructureParser();
        
        REQUIRE_THROWS_AS(p.parse("GG+CC", "(((+))"), sstruct::SecondaryStructureException);
        REQUIRE_THROWS_AS(p.parse("GG+CC", "()+))"), sstruct::SecondaryStructureException);
        REQUIRE_THROWS_AS(p.parse("GGG+CC", "(((+))"), sstruct::SecondaryStructureException);
        REQUIRE_THROWS_AS(p.parse("GG+CC", "(.+))"), sstruct::SecondaryStructureException);

    }
    
    SECTION("test parsing into graph") {
        auto p = sstruct::SecondaryStructureParser();
        auto g = p.parse("GGG+CC", ".((+))");
        REQUIRE(g->size() == 5);
        
        p.reset();
        
        SECTION("make sure that we are parsing the trailing unpaired residues") {
            g = p.parse("GG+CCAA", "((+))..");
            REQUIRE(g->size() == 5);
        }
        
    }
    
    SECTION("test parsing into single motifs") {
        auto p = sstruct::SecondaryStructureParser();
        
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
        auto p = sstruct::SecondaryStructureParser();
        auto motifs = p.parse_to_motifs("GGAGG+CAACCC", "((.((+)..)))");
        REQUIRE(motifs.size() == 3);
        REQUIRE(motifs[0]->mtype() == HELIX);
        REQUIRE(motifs[1]->mtype() == TWOWAY);
        
        p.reset();
        motifs = p.parse_to_motifs("GGAAGACAAGACAACC", "((..(.)..(.)..))");
        REQUIRE(motifs.size() == 4);
        REQUIRE(motifs[1]->mtype() == NWAY);

        p.reset();
        motifs = p.parse_to_motifs("GGAAGACAAGACAAGACAAGACCC",
                                   "((..(.)..(.)..(.)..(.)))");
        REQUIRE(motifs.size() == 6);
        REQUIRE(motifs[1]->mtype() == NWAY);


    }
    
    
}