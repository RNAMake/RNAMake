//
// Created by Joseph Yesselman on 10/22/17.
//


#include <iostream>
#include "../common.hpp"

#include <data_structure/graph.h>

TEST_CASE( "Test Graph Data Structure ", "[Graph]" ) {

    SECTION("test basic functionality") {

        auto g = data_structure::FixedEdgeDirectedGraph<int>();
        g.add_node(0, 3);
        g.add_node(1, 3, 0, data_structure::NodeIndexandEdge{0, 0});

        REQUIRE(g.get_num_nodes() == 2);
        REQUIRE(g.has_parent(0) == false);
        REQUIRE(g.has_parent(1) == true);

        g.setup_transversal(0);
        auto path = std::vector<Index>();
        for(auto const & n : g) {
            path.push_back(n->index());
        }
        REQUIRE(path.size() == 2);

        auto g2 = data_structure::FixedEdgeDirectedGraph<int>(g);
        g.remove_node(0);

        REQUIRE(g2.get_num_nodes() == 2);
        REQUIRE(g2.get_num_edges() == 1);

        g2.remove_edge(data_structure::NodeIndexandEdge{0, 0}, data_structure::NodeIndexandEdge{1, 0});

        REQUIRE(g2.get_num_edges() == 0);

        // test undirected graph
        auto g3 = data_structure::FixedEdgeUndirectedGraph<int>();
        g3.add_node(0, 3);
        g3.add_node(1, 3);
        g3.add_edge(data_structure::NodeIndexandEdge{0, 0}, data_structure::NodeIndexandEdge{1, 0});

    }

    SECTION("test new types") {

        struct X { int x; int y; };

        auto g1 = data_structure::FixedEdgeDirectedGraph<X>();
        g1.add_node(X{0, 0}, 3);

        REQUIRE(g1.get_node_data(0).x == 0);

        g1.setup_transversal(0);
        int i = 0;
        for(auto const & n : g1) {
            i ++;
        }
    }

    SECTION("test iteration forwarding") {
        struct GraphStruct {
        public:
            GraphStruct():
            graph_(data_structure::FixedEdgeDirectedGraph<int>()) {}

        public:
            typedef typename data_structure::FixedEdgeDirectedGraph<int>::const_iterator const_iterator;
            typedef typename data_structure::FixedEdgeDirectedGraph<int>::iterator iterator;

            iterator begin() { return graph_.begin(); }
            iterator end()   { return graph_.end(); }

            const_iterator begin() const { return graph_.begin(); }
            const_iterator end()   const { return graph_.end(); }

        private:
            data_structure::FixedEdgeDirectedGraph<int> graph_;

        };

        int i = 0;
        auto gs = GraphStruct();
        for(auto & n : gs) {
            i++;
        }


    }

}