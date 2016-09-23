//
//  graph_iterator_test.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>
#include <vector>

//RNAMake Headers
#include "data_structure/graph/graph.h"

void test_norm(
    std::vector<GraphStatic<int>> const & graphs,
    std::vector<int> & sum) {
    
    for(int i = 0; i < 10000; i++ ) {
        int count = 0;
        for(auto const & g : graphs) {
            for(auto const & n : g.nodes()) {
                auto n2 = n;
                count++;
            }
        }
        sum[i] = count;
    }
    
}

void test_iter(
    std::vector<GraphStatic<int>> const & graphs,
    std::vector<int> & sum) {
    
    for(int i = 0; i < 10000; i++ ) {
        int count = 0;
        for(auto const & g : graphs) {
            for(auto const & n : g) {
                auto n2 = n;
                count++;
            }
        }
        sum[i] = count;
    }
    
}


int main(int argc, const char * argv[]) {
    std::vector<GraphStatic<int>> graphs(100);
    std::vector<int> sum1(10000);
    std::vector<int> sum2(10000);

    for(int i = 0; i < 100; i++) {
        GraphStatic<int> g;
        for(int j = 0; j < i+100; j++) {
            g.add_data(j, -1, -1, -1, 2);
        }
        graphs[i] = g;
    }
    
    test_iter(graphs, sum2);
    test_norm(graphs, sum1);

    return 0;
}
