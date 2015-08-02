//
//  pdb_parser_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "pdb_parser_unittest.h"
#include "util/file_io.h"
#include "structure/structure.h"
#include "structure/structure_factory.h"


int
PDBParserUnittest::test_parse() {
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
    
    StructureFactory sf;
    StructureOP s2 = sf.get_structure(m_path);
    
    String path = unittest_resource_dir() + "/structure/test_str_to_structure.dat";
    Strings lines = get_lines_from_file(path);
    
    ResidueTypeSet rts;
    Structure s = str_to_structure(lines[0], rts);
    AtomOPs org_atoms = s.atoms();
    AtomOPs new_atoms = s2->atoms();
    
    auto org_res = s.residues();
    auto new_res = s2->residues();
    
    float dist;
    for(int i = 0; i < org_atoms.size(); i++) {
        if(org_atoms[i] == nullptr || new_atoms[i] == nullptr) { continue; }
        dist = org_atoms[i]->coords().distance(new_atoms[i]->coords());
        if(dist > 0.1) { return 0; }
    }
    
    return 1;
}

int
PDBParserUnittest::run() {
    
    if (test_parse() == 0)       { std::cout << "test_parse failed" << std::endl; }
    return 0;
    
    
}


void
PDBParserUnittest::run_all() {
    String name = "PDBParserUnittest";
    typedef int (PDBParserUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_parse"   ] = &PDBParserUnittest::test_parse;
  
    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
            if(result == 0) {
                std::cout << name << "::" << kv.first << " FAILED!" << std::endl;
            }
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
        }
        
    }
}


