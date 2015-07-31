//
//  tree_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/26/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "tree_unittest.h"
#include "data_structure/tree/tree_node.h"
#include "data_structure/tree/tree.h"


int
TreeUnittest::test_node() {
    
    auto n1 = std::make_shared<TreeNodeDynamic<int>>(1,0,0);
    auto n2 = std::make_shared<TreeNodeDynamic<int>>(1,0,0);
    n1->add_child(n2);
    n2->parent(n1);

    if(n2->parent()->data() != 1) {
        return 0;
    }
    
    return 1;
}

int
TreeUnittest::test_creation() {
    TreeDynamic<int> t;
    t.add_data(0);
    t.add_data(1);
    t.add_data(2, 0);
    
    try {
        t.get_node(99);
        throw std::runtime_error("unexpected_error");
    }
    catch(TreeException e) { }
    catch(...) { return 0; }
    
    if(t.get_node(2)->data() != 2) { return 0; }
    if(t.get_node(0)->children().size() != 2) { return 0; }
    
    TreeStatic<int> t1;
    t1.add_data(0, 2);
    t1.add_data(1, 2);
    t1.add_data(2, 2, 0);

    if(t1.get_node(2)->data() != 2) { return 0; }
    if(t1.get_node(0)->children().size() != 2) { return 0; }

    return 1;
}

int
TreeUnittest::test_add() {
    TreeDynamic<int> t;
    t.add_data(0);
    
    try {
        t.add_data(1, 1);
        throw std::runtime_error("unexpected_error");
    }
    catch(TreeException e) { }
    catch(...) { return 0; }
    
    TreeStatic<int> t1;
    t1.add_data(0, 2);
    t1.add_data(1, 2);
    t1.add_data(2, 2, 0);
    
    try {
        t1.add_data(1, 1, 0);
        throw std::runtime_error("unexpected_error");
    }
    catch(TreeException e) { }
    catch(...) { return 0; }
    
    return 1;
}

int
TreeUnittest::test_remove_node() {
    TreeDynamic<int> t;
    t.add_data(0);
    t.add_data(1);
    t.remove_node(0);
    if(t.size() != 1) { return 0; }
    if(t.get_node(1)->children().size() != 0) { return 0; }
    
    TreeStatic<int> t1;
    t1.add_data(0);
    t1.add_data(1);
    t1.remove_node(0);
    if(t1.size() != 1) { return 0; }
    int count = 0;
    for(auto const & c : t1.get_node(1)->children()) {
        if(c != nullptr) { count++; }
    }
    if(count != 0) {
        return 0;
    }
    
    return 1;
}

int
TreeUnittest::test_iter() {
    TreeDynamic<float> t;
    t.add_data(2.0f);
    t.add_data(3.0f);
    t.add_data(4.0f);

    int count = 0;
    for(auto const & n : t) {
        auto n1 = n;
        count++;
    }
    if(count != 3) { return 0; }
     
    return 1;
}


int
TreeUnittest::run() {
    if (test_node() == 0)          {  std::cout << "test_node failed" << std::endl; }
    if (test_creation() == 0)      {  std::cout << "test_creation failed" << std::endl; }
    if (test_add() == 0)           {  std::cout << "test_add failed" << std::endl; }
    if (test_remove_node() == 0)   {  std::cout << "test_remove_node failed" << std::endl; }
    if (test_iter() == 0)          {  std::cout << "test_iter failed" << std::endl; }
    return 1;
}

