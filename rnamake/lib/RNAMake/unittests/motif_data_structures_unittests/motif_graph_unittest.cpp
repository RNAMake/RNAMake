//
//  motif_graph_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <cxxabi.h>

#include "motif_graph_unittest.h"
#include "motif_data_structures/motif_graph.h"
#include "resources/resource_manager.h"
#include "build/build_motif_graph.h"
#include "build/build_motif_tree.h"

#include "secondary_structure_unittests/util.h"

namespace unittests {
namespace motif_structures {

void
MotifGraphUnittest::test_creation() {
    auto mg = MotifGraph();
    auto sterics = mg.get_bool_option("sterics");
    if(sterics != true) {
        throw UnittestException("was unable to retreieve correct option value");
    }
    
    mg.set_option_value("sterics", false);
    sterics = mg.get_bool_option("sterics");
    if(sterics != false) {
        throw UnittestException("was unable to retreieve correct option value");
    }
}

void
MotifGraphUnittest::test_add_motif() {
    auto mg = MotifGraph();
    auto m  = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    
    mg.add_motif(m);
    if(mg.size() != 1) {
        throw UnittestException("did not get the correct number of motifs");
    }
    
    mg.add_motif(m);
    if(mg.size() != 2) {
        throw UnittestException("did not get the correct number of motifs");
    }
    
    auto rna_struc = mg.get_structure();
    if(rna_struc->chains().size() != 2) {
        throw UnittestException("did not get the correct number of chains");
    }
    
    try {
        mg.add_motif("HELIX.IDEAL.600");
        throw "failed";
    }
    catch(MotifGraphException const & e) { }
    catch(...) {
        throw UnittestException("did not get expected error");
    }
}
    
void
MotifGraphUnittest::test_remove() {
    auto mg = MotifGraph();
    auto m  = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    mg.add_motif(m);
    mg.add_motif(m);
    
    mg.remove_motif(1);
    auto rna_struc = mg.get_structure();

    if(rna_struc->residues().size() != 4) {
        throw UnittestException("did not get the correct number of residues");
    }
    
    mg.remove_motif(0);
    
    rna_struc = mg.get_structure();
    
    if(rna_struc->residues().size() != 0) {
        throw UnittestException("did not get the correct number of residues");
    }
    
}
    
void
MotifGraphUnittest::test_remove_2() {
    auto m  = ResourceManager::getInstance().get_motif("HELIX.IDEAL.2");
    auto mg = MotifGraph();
    mg.add_motif(m);
    mg.add_motif(m);
    mg.add_motif(m);
    mg.remove_motif(1);
    auto rna_struc = mg.get_structure();
    if(rna_struc->chains().size() != 4) {
        throw UnittestException("did not get the correct number of chains");
    }
    
    int count = 0;
    for(auto const & n : mg) {
        count++;
        auto n2 = n;
    }

    if(count != 2) {
        throw UnittestException("did not iterate through all nodes");
    }
    
}

void
MotifGraphUnittest::test_copy() {
    auto builder = BuildMotifGraph();
    auto mg = builder.build(3);
    
    auto mg_copy = std::make_shared<MotifGraph>(*mg);
    mg_copy->remove_motif(1);
    
    if(mg_copy->size() != 2) {
        throw UnittestException("did not get the correct number of motifs");
    }
    
    auto rna_struc = mg_copy->get_structure();
    if(rna_struc->chains().size() != 4) {
        throw UnittestException("did not get the correct number of chains");
    }
}
    
void
MotifGraphUnittest::test_replace_ideal_helices() {
    auto mg = MotifGraph();
    mg.add_motif("HELIX.IDEAL.6");
    mg.replace_ideal_helices();
    
    if(mg.size() != 7) {
        throw UnittestException("did not get the correct number of motifs");
    }

    auto builder = BuildMotifGraph();
    auto mg2 = builder.build(3);
    mg2->replace_ideal_helices();
}
    
void
MotifGraphUnittest::test_replace_ideal_helices_2() {
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/mini_ttr/resources/";
    ResourceManager::getInstance().add_motif(base_path+"GAAA_tetraloop");
    
    auto mg = MotifGraph();
    mg.add_motif("GAAA_tetraloop", "A229-A245");
    mg.add_motif("HELIX.IDEAL.6", -1, "A149-A154");
    mg.replace_ideal_helices();

}
    
void
MotifGraphUnittest::test_secondary_structure() {
    auto builder = BuildMotifGraph();
    auto mg = builder.build(3);
    auto ss = mg->secondary_structure();
    if(ss->motifs().size() != 3) {
        throw UnittestException("did not get the correct number of motifs");
    }

    mg->replace_ideal_helices();
    ss = mg->secondary_structure();
    if(mg->size() != ss->motifs().size()) {
        throw UnittestException("did not get the correct number of motifs");
    }
    
}
    
void
MotifGraphUnittest::test_replace_sequence() {
    auto builder = BuildMotifGraph();
    auto mg = builder.build(3);
    mg->replace_ideal_helices();
    auto ss = mg->designable_secondary_structure();
    unittests::sstruct_unittests::fill_basepairs_in_ss(ss);
    mg->replace_helical_sequence(ss);
}
    

//no leak from just copying
void
MotifGraphUnittest::test_memory() {
    auto builder = BuildMotifTree();
    auto mg = builder.build(20);
    int count = 0;
    for(int i = 0; i < 10000; i++) {
        auto mg2 = std::make_shared<MotifTree>();
        for(auto const & n : *mg) {
            mg2->add_motif(n->data());
        }
        count += mg2->size();
    }
}

//no leak
void
MotifGraphUnittest::test_memory_2() {
    
    auto builder = BuildMotifTree();
    auto mlib = std::make_shared<MotifSqliteLibrary>("twoway");
    for(int i = 0; i < 10000; i++) {
        auto mg = std::make_shared<MotifGraph>();
        for(int j = 0; j < 20; j++) {
            auto m = mlib->get_random();
        }
    }
}
    
void
MotifGraphUnittest::test_memory_3() {
    auto builder = BuildMotifGraph();
    for(int i = 0; i < 10000; i++) {
        auto mg = builder.build(20);
        //std::cout << mg->size() << std::endl;
    }
}
    
int
MotifGraphUnittest::run() {
    test_remove();

    /*test_creation();
    test_add_motif();
    test_remove_2();
    test_copy();
    test_replace_ideal_helices();
    test_replace_ideal_helices_2();
    test_secondary_structure();
    test_replace_sequence();*/
    //test_memory_2();
    return 0;
}
    
int
MotifGraphUnittest::run_all() {
    String name = "MotifGraphUnittest";
    typedef void (MotifGraphUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"      ] = &MotifGraphUnittest::test_creation;
    func_map["test_add_motif"     ] = &MotifGraphUnittest::test_add_motif;
    func_map["test_remove"        ] = &MotifGraphUnittest::test_remove;
    func_map["test_remove_2"      ] = &MotifGraphUnittest::test_remove_2;
    func_map["test_copy"          ] = &MotifGraphUnittest::test_copy;
    func_map["test_replace_ideal_helices"  ] = &MotifGraphUnittest::test_replace_ideal_helices;
    func_map["test_replace_ideal_helices_2"] = &MotifGraphUnittest::test_replace_ideal_helices_2;
    func_map["test_secondary_structure"    ] = &MotifGraphUnittest::test_secondary_structure;
    func_map["test_replace_sequence"       ] = &MotifGraphUnittest::test_replace_sequence;


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



















