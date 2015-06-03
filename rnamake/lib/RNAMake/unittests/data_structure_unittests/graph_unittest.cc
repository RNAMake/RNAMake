//
//  graph_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "graph_unittest.h"
#include "graph.h"

int
GraphUnittest::test_creation() {
    Graph<int> g;
    g.add_data(1);

    return 1;
}



int
GraphUnittest::run() {
    if (test_creation() == 0)            { std::cout << "test_creation failed" << std::endl;  }
    return 0;
}