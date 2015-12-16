//
//  secondary_structure_parser_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_parser_unittest__
#define __RNAMake__secondary_structure_parser_unittest__

#include <stdio.h>


#include "unittest.h"

namespace unittests {
    
class SecondaryStructureParserUnittest : public Unittest {
public:
    SecondaryStructureParserUnittest() {}
    
    ~SecondaryStructureParserUnittest() {}
    
public:
    
    void
    test_creation();
    
    void
    test_parse();
    
    void
    test_parse_to_motifs();
    
    void
    test_tecto();
    
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};
    
void
debug_ss_parser_graph(
    String const &,
    String const &);
    
}




#endif /* defined(__RNAMake__secondary_structure_parser_unittest__) */
