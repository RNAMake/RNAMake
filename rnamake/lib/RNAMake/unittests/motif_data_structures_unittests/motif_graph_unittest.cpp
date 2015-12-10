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

namespace unittests {
namespace motif_structures {

void
MotifGraphUnittest::test_creation() {
    auto mg = MotifGraph();
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
    catch(MotifMergerException const & e) { }
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
MotifGraphUnittest::test_secondary_structure() {
    auto builder = BuildMotifGraph();
    auto mg = builder.build(3);
    auto ss = mg->secondary_structure();
    if(ss->motifs().size() != 3) {
        throw UnittestException("did not get the correct number of motifs");
    }

    mg->replace_ideal_helices();
    ss = mg->secondary_structure();
    if(mg->size() == ss->motifs().size()) {
        throw UnittestException("did not get the correct number of motifs");
    }
    
}
    
    
int
MotifGraphUnittest::run() {
    test_creation();
    test_add_motif();
    test_remove();
    test_remove_2();
    test_copy();
    test_replace_ideal_helices();
    test_secondary_structure();
    return 0;
}
    
}
}



















