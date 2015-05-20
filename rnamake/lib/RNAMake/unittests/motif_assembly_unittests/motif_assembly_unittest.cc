//
//  motif_assembly_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_assembly_unittest.h"
#include "resources/motif_library.h"
#include "resources/library_manager.h"
#include "motif_assembly/motif_assembly.h"

int
MotifAssemblyUnittest::test_add_motif() {
    MotifLibrary mlib(HELIX);
    MotifAssembly ma;
    ma.add_motif(mlib.get_motif("HELIX.IDEAL"));

    if(ma.nodes().size() != 2) { return 0; }
    
    MotifAssembly ma2;
    ma2.add_motif("HELIX.IDEAL");
    
    if(ma2.nodes().size() != 2) { return 0; }

    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/fmn_min";
    LibraryManager::getInstance().add_motif(m_path);
    
    MotifAssembly ma3;
    ma3.add_motif("fmn_min");
    
    return 1;
}

int
MotifAssemblyUnittest::run() {
    if (test_add_motif() == 0)              {  std::cout << "test_add_motif failed" << std::endl; }
    //if (test_motif_end_indentity() == 0)    {  std::cout << "test_motif_end_indentity failed" << std::endl; }

    return 1;
}