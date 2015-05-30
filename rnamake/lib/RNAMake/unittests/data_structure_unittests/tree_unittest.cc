//
//  tree_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/26/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "tree_unittest.h"
#include "data_structure/tree/node.h"
#include "data_structure/tree/tree.h"


int
TreeUnittest::test_node() {
    
    auto n1 = new Node<int>(1,0,0);
    auto n2 = new Node<int>(1,0,0);
    n1->add_child(n2);
    n2->set_parent(n1);

    if(n2->parent()->data() != 1) {
        return 0;
    }
    
    delete n1;
    delete n2;
    
    return 1;
}

int
TreeUnittest::test_creation() {
    Tree <int> t;
    t.add_data(0);
    t.add_data(1);
    t.add_data(2, 0, 0);
    
    if (t.get_data(2) != 2) {
        return 0;
    }
    
    if(t.get_node(0)->children().size() != 2) { return 0; }
    
    Tree<int> t2(NodeTypeStatic);
    t2.add_data(10, 1);
    t2.add_data(11, 1, 0);
    try {
        t2.add_data(12, 1, 0);
        std::cout << "did not catch error properly\n";
        exit(EXIT_FAILURE);
    } catch(...) {}

    return 1;
}

int
TreeUnittest::test_remove_node() {
    Tree <int> t;
    t.add_data(0);
    t.add_data(1);
    t.add_data(2, 0, 0);
    t.add_data(3, 0, 1);
    t.remove_node(0);
    if(t.size() != 0) { return 0; }
    return 1;
}

int
TreeUnittest::test_get_index() {
    Tree<float> t;
    t.add_data(2.0f);
    t.add_data(3.0f);
    t.add_data(4.0f);
    
    int index = t.get_index(3.0f);
    if(index != 1) { return 0; }
    
    return 1;
}

int
TreeUnittest::run() {
    if (test_node() == 0)          {  std::cout << "test_node failed" << std::endl; }
    if (test_creation() == 0)      {  std::cout << "test_creation failed" << std::endl; }
    if (test_remove_node() == 0)   {  std::cout << "test_remove_node failed" << std::endl; }
    if (test_get_index() == 0)     {  std::cout << "test_get_index failed" << std::endl; }

    return 1;
}