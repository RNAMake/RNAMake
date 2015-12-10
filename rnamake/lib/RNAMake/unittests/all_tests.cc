//
//  all_tests.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>
#include <vector>

//RNAMake Headers
#include "base_unittests/option_unittest.h"

#include "data_structure_unittests/graph_unittest.h"



int main(int argc, const char * argv[]) {
    std::vector<Unittest*> units;
    units.push_back(new unittests::OptionUnittest());
    units.push_back(new unittests::CL_OptionUnittest());

    //units.push_back(new GraphUnittest());
    /*units.push_back(new UuidUnittest());
    units.push_back(new X3dnaUnittest());
    units.push_back(new AtomUnittest());
    units.push_back(new ResidueTypeSetUnittest());
    units.push_back(new ResidueUnittest());
    units.push_back(new ChainUnittest());
    units.push_back(new StructureUnittest());
    units.push_back(new BasepairUnittest());
    units.push_back(new PDBParserUnittest());
    units.push_back(new MotifUnittest());
    units.push_back(new MotifTreeUnittest());
    units.push_back(new MotifTreeMergerUnittest());
    units.push_back(new MotifScorerUnittest());
    units.push_back(new MotifTreeStateLibraryUnittest());
    units.push_back(new MotifLibraryUnittest());
    units.push_back(new LibraryManagerUnittest());*/
    
    int tests_run = 0;
    int test_failed = 0;
    for(auto const & test : units) {
        if(test == nullptr) { continue; }
        test_failed += test->run_all();
        test_run    += test->size();
    }
    std::cout "RAN: " << test_run << "  FAILED: " << test_failed << std::endl;

    
    
    return 0;
}












