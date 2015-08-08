//
//  motif_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "unittest.h"
#include "motif_unittest.h"
#include "motif/motif_factory.h"
#include "motif/motif_to_secondary_structure.h"
#include "util/file_io.h"
#include "util/settings.h"
#include "math/numerical.h"

MotifUnittest::MotifUnittest() {
    
    MotifFactory mf;
    String path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6";
    p4p6_ = mf.motif_from_file(path);

    String path1 = base_dir() + "/rnamake/resources/motifs/base.motif";
    base_ = file_to_motif(path1);
    
}

int
MotifUnittest::test_copy() {
    Motif mcopy = p4p6_->copy();
    Point p(50,0,0);
    mcopy.move(p);
    
    Point org_center = center(p4p6_->atoms());
    Point new_center = center(mcopy.atoms());
    float dist = (50) - org_center.distance(new_center);
    if(dist > 0.001) { return 0; }
    
    return 1;
}

int
MotifUnittest::test_to_str() {
    ResidueTypeSet rts;
    Motif m2 = Motif(p4p6_->to_str(), rts);
    
    AtomOPs new_atoms = m2.atoms();
    AtomOPs org_atoms = p4p6_->atoms();
    
    float dist;
    for(int i = 0; i < org_atoms.size(); i++) {
        if(org_atoms[i].get() == NULL) { continue; }
        dist = org_atoms[i]->coords().distance(new_atoms[i]->coords());
    }
    return 1;
}

int
MotifUnittest::test_secondary_structure() {
    String org = "....((((((...(((.((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....)).)))..))))))(....((((..(((((((((.....)))))))).)))))......)";
    String new_ss = p4p6_->dot_bracket();
    for(int i = 0; i < org.size(); i++) {
        if(org[i] != new_ss[i]) { return 0; }
    }
    return 1;
}

int
MotifUnittest::test_get_basepair_by_name() {
    String path = base_dir() + "/rnamake/resources/motifs/base.motif";
    auto m = file_to_motif(path);
    auto bp = m->get_basepair("A2-A9")[0];

    return 1;
}

int
MotifUnittest::test_align() {
    auto m = std::make_shared<Motif>(base_->copy());
    auto m1 = std::make_shared<Motif>(base_->copy());
    
    align_motif(m->ends()[1]->state(), m1->ends()[0], m1);
    
    if(!are_xyzVector_equal(m->ends()[1]->d(), m1->ends()[0]->d())) {
        return 0;
    }
        
    return 1;
}


int
MotifUnittest::run() {
    if (test_copy() == 0)                 { std::cout << "test_copy failed" << std::endl; }
    if (test_to_str() == 0)               { std::cout << "test_to_str failed" << std::endl; }
    if (test_secondary_structure() == 0)  { std::cout << "test_secondary_structure failed" << std::endl; }
    //if (test_get_basepair_by_name() == 0) { std::cout << "test_get_basepair_by_name failed" << std::endl; }
    if (test_align() == 0)                { std::cout << "test_align failed" << std::endl; }

    return 0;
}

void
MotifUnittest::run_all() {
    /*String name = "MotifUnittest";
    typedef int (MotifUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_copy"   ] = &MotifUnittest::test_copy;
    func_map["test_to_str" ] = &MotifUnittest::test_to_str;
    func_map["test_secondary_structure"] = &MotifUnittest::test_secondary_structure;
    func_map["test_creation_from_dir"       ] = &MotifUnittest::test_creation_from_dir;
    func_map["test_get_basepair_by_name"       ] = &MotifUnittest::test_get_basepair_by_name;

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
























