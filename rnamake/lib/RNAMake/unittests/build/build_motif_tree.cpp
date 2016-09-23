//
//  build_motif_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "build/build_motif_tree.h"

MotifTreeOP
BuildMotifTree::build(
    int size) {
    
    auto mt = std::make_shared<MotifTree>();
    int i = 0, pos = 0, count = 0;
    while ( mt->size() < size) {
        if(pos == libs_.size()) { pos = 0; }
        
        auto m = libs_[pos]->get_random();
        i = mt->add_motif(m);
        count += 1;
        if (count > 100) { break; }
        if (i != -1) { pos += 1; }
        else {
            std::cout << "fail: " << m->name() << std::endl;
            return mt;
        }
    }
    
    return mt;
}


std::unique_ptr<MotifTree>
BuildMotifTree::build_2(
    int size) {
    
    auto mt = std::make_unique<MotifTree>();
    //mt->option("sterics", 0);
    int i = 0, pos = 0, count = 0;
    while ( mt->size() < size) {
        if(pos == libs_.size()) { pos = 0; }
        
        auto m = libs_[pos]->get_random();
        //i = mt->add_motif(m);
        count += 1;
        if (count > 100) { break; }
        if (i != -1) { pos += 1; }
        else {
            std::cout << "fail: " << m->name() << std::endl;
            return mt;
        }
    }
    
    return mt;
}

