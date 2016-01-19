//
//  graph_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>
#include "graph_unittest.h"
#include "data_structure/graph/graph.h"


int
GraphUnittest::test_nodes() {
    auto n1 = std::make_shared<GraphNodeDynamic<int>>(0, 0, 0);
    auto n2 = std::make_shared<GraphNodeDynamic<int>>(1, 1, 0);
    auto c1 = std::make_shared<GraphConnection<int>>(n1, n2, 0, 0);
    n1->add_connection(c1);
    n2->add_connection(c1);
    
    auto n3 = std::make_shared<GraphNodeStatic<int>>(0, 0, 0, 2);
    auto n4 = std::make_shared<GraphNodeStatic<int>>(1, 1, 0, 2);
    auto c2 = std::make_shared<GraphConnection<int>>(n1, n2, 1, 1);
    n3->add_connection(c2, 1);
    n4->add_connection(c2, 1);

    try {
        n3->add_connection(c2, 2);
        throw UnittestException("did not catch GraphException");
    } catch (GraphException & e) {}
    
    Ints avail_pos = n3->available_children_pos();
    
    if(avail_pos.size() != 1 || avail_pos[0] != 0) {
        return 0;
    }
    
    return 1;
}

int
GraphUnittest::test_creation() {
    
    GraphDynamic<int> g;
    g.add_data(0);
    g.add_data(1);
    g.add_data(2, 0);
    g.add_data(3, 0);
    
    //try to get a node that doesnt exist
    try { g.get_node(10); throw std::runtime_error("failed"); }
    catch(GraphException e) {}
    catch(...) { return 0; }
    
    
    GraphStatic<int> g1;
    g1.add_data(0, -1, -1, -1, 2);
    g1.add_data(1, -1, -1, -1, 2);
    
    return 1;
}

int
GraphUnittest::test_add() {
    GraphDynamic<int> g;
    g.add_data(0);
    
    //catch improper parent index
    try { g.add_data(1, 1); throw std::runtime_error("failed"); }
    catch(GraphException e) {}
    catch(...) { return 0; }
    
    GraphStatic<int> g1;
    g1.add_data(0, -1, -1, 0, 1);
    
    //catch improper parent index
    try { g1.add_data(1, 1); throw std::runtime_error("failed"); }
    catch(GraphException e) {}
    catch(...) { return 0; }

    g1.add_data(1, 0, 0, 0, 2);
    
    //catch improper connection index, cannot add to node 0 already at max connections
    try { g1.add_data(2, 0); throw std::runtime_error("failed");  }
    catch(GraphException e) {}
    catch(...) { return 0; }
    
    //catch incorrect end index,
    try { g1.add_data(2, 1, 0, 0, 1); throw std::runtime_error("failed"); }
    catch(GraphException e) {}
    catch(...) { return 0; }
    
    g1.add_data(2, 1, 1, 0, 1);
    
    return 1;
    
}

int
GraphUnittest::test_connect() {
    GraphStatic<int> g;
    g.add_data(0, -1, -1, -1, 2);
    g.add_data(1,  0, 0, 0, 2);
    g.add_data(2,  1, 1, 0, 2);
    g.connect(0, 2, 1, 1);
    return 1;
}

int
GraphUnittest::test_remove() {
    GraphStatic<int> g;
    g.add_data(0, -1, -1, -1, 2);
    g.remove_node(0);
    
    g.add_data(0, -1, -1, -1, 2);
    g.add_data(1,  1, 0, 0, 2);
    g.add_data(2,  2, 1, 0, 2);
    g.remove_node(3);
    
    if(g.size() != 2) {
        return 0;
    }
    
    GraphStatic<int> g1;
    g1.add_data(0, -1, -1, -1, 2);
    g1.add_data(1,  0, 0, 0, 2);
    g1.add_data(2,  1, 1, 0, 2);
    g1.remove_node(0);
    
    if(g.size() != 2) {
        return 0;
    }
    
    return 1;
}

int
GraphUnittest::test_iteration() {
    GraphStatic<int> g;
    g.add_data(0, -1, -1, -1, 2);
    g.add_data(1,  0, 0, 0, 2);
    g.add_data(2,  1, 1, 0, 2);
    int count = 0;
    for(auto const & n : g) { auto n2 = n; count++; }
    if(count != 3) {
        return 0;
    }
    count = 0;
    g.remove_node(0);
    for(auto const & n : g) { auto n2 = n; count++; }
    if(count != 2) {
        return 0;
    }
    
    return 1;
}


int
GraphUnittest::run() {
    if (test_nodes() == 0)           { std::cout << "test_nodes failed" << std::endl;  }
    if (test_creation() == 0)        { std::cout << "test_creation failed" << std::endl;  }
    if (test_add() == 0)             { std::cout << "test_add failed" << std::endl;  }
    if (test_connect() == 0)         { std::cout << "test_connect failed" << std::endl;  }
    if (test_remove() == 0)          { std::cout << "test_remove failed" << std::endl;  }
    if (test_iteration() == 0)       { std::cout << "test_iteration failed" << std::endl;  }

    return 0;
}

int
GraphUnittest::run_all() {
    String name = "GraphUnittest";
    typedef int (GraphUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_nodes"      ] = &GraphUnittest::test_nodes;
    func_map["test_creation"   ] = &GraphUnittest::test_creation;
    func_map["test_add"        ] = &GraphUnittest::test_add;
    func_map["test_connect"    ] = &GraphUnittest::test_connect;
    func_map["test_remove"     ] = &GraphUnittest::test_remove;
    func_map["test_iteration"  ] = &GraphUnittest::test_iteration;

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
        std::cout << ".";
        
    }
    
    return 0;
}