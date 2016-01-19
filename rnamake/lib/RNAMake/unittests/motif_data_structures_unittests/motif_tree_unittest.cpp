//
//  motif_tree_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_unittest.h"
#include "motif_data_structures/motif_tree.h"
#include "resources/resource_manager.h"

namespace unittests {
namespace motif_structures {

void
MotifTreeUnittest::test_creation() {
    auto mt = MotifTree();
    auto sterics = mt.get_bool_option("sterics");
    if(sterics != true) {
        throw UnittestException("was unable to retreieve correct option value");
    }
    
    mt.set_option_value("sterics", false);
    sterics = mt.get_bool_option("sterics");
    if(sterics != false) {
        throw UnittestException("was unable to retreieve correct option value");
    }

}
    
void
MotifTreeUnittest::test_add_motif() {
    auto mt = MotifTree();
    auto m  = ResourceManager::getInstance().get_motif("HELIX.IDEAL.2");
    mt.add_motif(m);
    if(mt.size() != 1) {
        throw UnittestException("did not add motif properly to tree");
    }
    
    mt.add_motif(m);
    if(mt.size() != 2) {
        throw UnittestException("did not add motif properly to tree");
    }
    
    auto rna_struc = mt.get_structure();
    if(rna_struc->residues().size() != 14) {
        throw UnittestException("merger did not result in the right number of residues");
    }
}

    
int
MotifTreeUnittest::run() {
    test_creation();
    test_add_motif();
    return 0;
}
    
    
int
MotifTreeUnittest::run_all() {
    String name = "MotifTreeUnittest";
    typedef void (MotifTreeUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"      ] = &MotifTreeUnittest::test_creation;
    func_map["test_add_motif"     ] = &MotifTreeUnittest::test_add_motif;
    
    
    int failed = 0;
    for(auto const & kv : func_map) {
        try {
            (this->*kv.second)();
        }
        catch(std::exception const & e) {
            std::cout << name << "::" << kv.first << " returned ERROR! : " << e.what() << std::endl;
            failed++;
        }
    }
    return failed;
}
    
    
}
}