//
//  graph_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "graph_unittest.h"
//#include "data_structure/graph/graph.h"
#include "data_structure/graph/graph_node.h"


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
    n3->add_connection(c2, 0);
    n4->add_connection(c2, 0);

    /*try {

    } catch (GraphException & e) { }
    */
    
    
    
    return 1;
    
}

int
GraphUnittest::test_creation() {
    /*Graph<int> g;
    g.add_data(0);
    g.add_data(1);
    g.add_data(2, 0, 0);
    
    auto n = g.get_node(2);
    
    if(n->parent()->data() != 0) { return 0; }
    
    Graph<int> g2(GraphNodeType::GraphNodeTypeStatic);
    g2.add_data(0, 2);
    */
    
    
    return 1;
}



int
GraphUnittest::run() {
    if (test_nodes() == 0)               { std::cout << "test_nodes failed" << std::endl;  }
    if (test_creation() == 0)            { std::cout << "test_creation failed" << std::endl;  }
    return 0;
}