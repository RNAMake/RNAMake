//
//  motif_tree_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_unittest.h"
#include "motif/motif_tree.h"
#include "resources/resource_manager.h"

int
MotifTreeUnittest::test_creation() {
    MotifTree mt;
    return 1;
}

int
MotifTreeUnittest::test_add_motif() {
    MotifTree mt;
    MotifOP m = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    for (int i = 0; i < 20; i ++) { mt.add_motif(m); }
    std::cout << mt.size() << std::endl;
    mt.write_pdbs();
    return 1;
    
    try {
        mt.add_motif(m, 0);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(char const * e) { }
    
    try {
        mt.add_motif(m, -1, 10);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(char const * e) { }
    
    
    try {
        mt.add_motif(m, -1, -1);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(char const * e) { }
    

    return 1;
}

/*
int
MotifTreeUnittest::test_motif_tree_to_str() {
    String path = unittest_resource_dir() + "/motif_tree/test_motif_tree_to_str.dat";
    Strings lines = get_lines_from_file(path);
    ResidueTypeSet rts;
    MotifLibrary mlib ( TWOWAY );
    int i = 0;
    for (auto const & line : lines) {
        if(line.length() < 5) { break; }
        MotifTree mt ( line, rts );
        MotifTree mt2;
        int j = -1;
        for (auto const & n : mt.nodes()) {
            j++;
            if (j == 0) { continue; }
            MotifOP m = mlib.get_motif(n->motif()->name());
            MotifTreeNodeOP node = mt2.add_motif(m, nullptr, 0);
            if(node == nullptr) {
                std::cout << m->name() << std::endl;
            }
        }
        Point org = mt.nodes().back()->motif()->ends()[1]->d();
        Point np = mt2.nodes().back()->motif()->ends()[1]->d();
        float dist = org.distance(np);
        if (dist > 0.1) {
            return 0;
        }
        i++;
        
    }
    return 1;
}

int
MotifTreeUnittest::test_remove_node() {
    MotifLibrary mlib (HELIX);
    MotifTree mt;
    MotifOP m = mlib.get_motif("HELIX.IDEAL");
    mt.add_motif(m);
    mt.remove_node(mt.last_node());
    
    return 1;
}

int
MotifTreeUnittest::test_remove_node_level() {
    MotifLibrary mlib (HELIX);
    MotifTree mt;
    MotifOP m = mlib.get_motif("HELIX.IDEAL");
    for(int i = 0; i < 10; i++) {
        mt.add_motif(m);
    }
    mt.remove_node_level();
    if(mt.nodes().size() != 1) { return 0; }
    return 1;
}

int
MotifTreeUnittest::test_options() {
    MotifLibrary mlib (HELIX);
    MotifTree mt;
    MotifOP m = mlib.get_motif("HELIX.IDEAL");
    
    if(mt.option<float>("clash_radius") != 2.8f) { return 0; }
    
    mt.option("clash_radius", 3.0f);
    
    if(mt.option<float>("clash_radius") != 3.0f) { return 0; }

    
    return 1;
}
*/
 
int
MotifTreeUnittest::run() {
    if (test_creation() == 0)          { std::cout << "test_creation failed" << std::endl;  }
    if (test_add_motif() == 0)         { std::cout << "test_add_motif failed" << std::endl; }
    //if (test_motif_tree_to_str() == 0) { std::cout << "test_motif_tree_to_str failed" << std::endl; }
    //if (test_remove_node() == 0)       { std::cout << "test_remove_node failed" << std::endl; }
    //if (test_remove_node_level() == 0) { std::cout << "test_remove_node_level failed" << std::endl; }
    //if (test_options() == 0) { std::cout << "test_options failed" << std::endl; }
    
    return 0;
}


void
MotifTreeUnittest::run_all() {
    /*String name = "MotifTreeUnittest";
    typedef int (MotifTreeUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"          ] = &MotifTreeUnittest::test_creation;
    func_map["test_add_motif"         ] = &MotifTreeUnittest::test_add_motif;
    func_map["test_motif_tree_to_str" ] = &MotifTreeUnittest::test_motif_tree_to_str;
    func_map["test_remove_node"       ] = &MotifTreeUnittest::test_remove_node;
    func_map["test_remove_node_level" ] = &MotifTreeUnittest::test_remove_node_level;
    func_map["test_options"           ] = &MotifTreeUnittest::test_options;
    
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
        
    }*/
}

