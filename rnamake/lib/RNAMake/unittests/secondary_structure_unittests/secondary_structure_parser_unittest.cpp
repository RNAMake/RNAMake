//
//  secondary_structure_parser_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_parser_unittest.h"
#include "secondary_structure/secondary_structure_parser.h"

namespace unittests {
    
void
SecondaryStructureParserUnittest::test_creation() {
    sstruct::SecondaryStructureParser parser;
}

void
SecondaryStructureParserUnittest::test_parse() {
    sstruct::SecondaryStructureParser parser;
    auto g = parser.parse("GG+CC", "((+))");
    if(g->size() != 4) {
        throw UnittestException("did not get the right number of nodes");
    }
    
    /*for(auto const & n : *g) {
        std::cout << n->index() << " ";
        for(auto const & c : n->connections()) {
            if(c == nullptr) { std::cout << "N "; }
            else { std::cout << c->partner(n->index())->index() << " "; }
        }
        std::cout << "\n";
    }*/
    
    auto g1 = parser.parse("GAAG+CAC", "(..(+).)");
    if(g1->size() != 6) {
        throw UnittestException("did not get the right number of nodes");
    }


}
    
void
SecondaryStructureParserUnittest::test_parse_to_motifs() {
    sstruct::SecondaryStructureParser parser;
    auto motifs = parser.parse_to_motifs("GG+CC", "((+))");
    if(motifs.size() != 1) {
        throw UnittestException("did not get the correct number of motifs");
    }
}

int
SecondaryStructureParserUnittest::run() {
    test_creation();
    test_parse();
    test_parse_to_motifs();
    return 0;
}

    
}