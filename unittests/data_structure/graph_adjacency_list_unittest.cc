//
// Created by Joseph Yesselman on 10/22/17.
//


#include <iostream>
#include "../common.hpp"

#include <data_structure/graph_adjacency_list.h>

TEST_CASE( "Test Graph Data Structure ", "[Graph]" ) {

    SECTION("testing adding nodes and edges") {

        auto adj_list = data_structure::FixedEdged_AL<int>();
        adj_list.add_node(0, 3);
        adj_list.add_node(1, 3);

        REQUIRE(adj_list.get_num_nodes() == 2);
        REQUIRE(adj_list.get_num_edges() == 0);

        adj_list.add_edge(data_structure::NodeIndexandEdge{0, 0},
                          data_structure::NodeIndexandEdge{1, 0});

        REQUIRE(adj_list.get_num_edges() == 1);

        adj_list.remove_node(0);

        REQUIRE(adj_list.get_num_nodes() == 1);
        REQUIRE(adj_list.get_num_edges() == 0);

        adj_list = data_structure::FixedEdged_AL<int>();
        adj_list.add_node(0, 3);
        adj_list.add_node(1, 3);

        adj_list.add_edge(data_structure::NodeIndexandEdge{0, 0},
                          data_structure::NodeIndexandEdge{1, 0});

        adj_list.remove_edge(data_structure::NodeIndexandEdge{0, 0},
                             data_structure::NodeIndexandEdge{1, 0});

        REQUIRE(adj_list.get_num_edges() == 0);
    }

    SECTION("test copying") {

        auto adj_list = data_structure::FixedEdged_AL<int>();
        adj_list.add_node(0, 3);
        adj_list.add_node(1, 3);
        adj_list.add_edge(data_structure::NodeIndexandEdge{0, 0},
                          data_structure::NodeIndexandEdge{1, 0});

        auto adj_list_2 = adj_list;

        adj_list.remove_node(0);

        REQUIRE(adj_list_2.get_num_edges() == 1);
        REQUIRE(adj_list_2.get_num_nodes() == 2);

    }

    SECTION("test copying 2") {
        auto adj_list = data_structure::FixedEdged_AL<int>();
        adj_list.add_node(0, 3);
        adj_list.add_node(1, 3);
        adj_list.add_edge(data_structure::NodeIndexandEdge{0, 0},
                          data_structure::NodeIndexandEdge{1, 0});

        auto adj_list_2 = data_structure::FixedEdged_AL<int>(adj_list);

        adj_list.remove_node(0);

        REQUIRE(adj_list_2.get_num_edges() == 1);
        REQUIRE(adj_list_2.get_num_nodes() == 2);

    }

    SECTION("test directed") {
        auto adj_list = data_structure::FixedEdged_DAL<int>();
        adj_list.add_node(0, 3);
        adj_list.add_node(1, 3, 0, data_structure::NodeIndexandEdge{0, 0});

        REQUIRE(adj_list.get_num_nodes() == 2);

        REQUIRE(adj_list.has_parent(0) == false);
        REQUIRE(adj_list.has_parent(1) == true);

    }

    SECTION("test directed copy") {
        auto adj_list = data_structure::FixedEdged_DAL<int>();
        adj_list.add_node(0, 3);
        adj_list.add_node(1, 3, 0, data_structure::NodeIndexandEdge{0, 0});

        auto adj_list_2 =  data_structure::FixedEdged_DAL<int>(adj_list);

        adj_list.remove_node(0);

        REQUIRE(adj_list_2.get_num_nodes() == 2);

        REQUIRE(adj_list_2.has_parent(0) == false);
        REQUIRE(adj_list_2.has_parent(1) == true);

        adj_list_2.add_node(2, 3, 0, data_structure::NodeIndexandEdge{1, 1});

        REQUIRE(adj_list_2.get_num_nodes() == 3);

    }

}