//
//  pdb_parser_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__pdb_parser_unittest__
#define __RNAMake__pdb_parser_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"
#include "structure/pdb_parser.h"

class PDBParserUnittest : public Unittest {
public:
    
    PDBParserUnittest() {}
    
    ~PDBParserUnittest() {}
    
public:
    
    int
    test_parse();
    
public:
    
    int
    run();
    
private:
    
    
};

#endif /* defined(__RNAMake__pdb_parser_unittest__) */
