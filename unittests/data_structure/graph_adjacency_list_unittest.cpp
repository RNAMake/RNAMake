//
// Created by Joseph Yesselman on 10/22/17.
//

#include "../common.hpp"
#include <iostream>

#include <data_structure/graph/graph_adjacency_list.h>

TEST_CASE("Test Graph Data Structure ") {
  using namespace data_structure::graph;
  SUBCASE("test basic adjacency list") {
    SUBCASE("test add") {
      int i = 0;
      AdjacencyList<int, FixedEdges> adj_list;
      adj_list.add_node(i, 1);
      CHECK(adj_list.get_num_connections() == 0);
      CHECK(adj_list.get_num_nodes() == 1);
      CHECK(adj_list.get_node_data(0) == 0);
    }
    SUBCASE("test change value of data") {
      int i = 0, j = 1;
      AdjacencyList<int, FixedEdges> adj_list;
      adj_list.add_node(i, 1);
      adj_list.set_node_data(0, j);
      CHECK(adj_list.get_node_data(0) == 1);
    }
    SUBCASE("test add connection") {
      int i = 0, j = 1;
      AdjacencyList<int, FixedEdges> adj_list;
      adj_list.add_node(i, 1);
      adj_list.add_node(j, 1);
      adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
      CHECK(adj_list.get_num_nodes() == 2);
      CHECK(adj_list.get_num_connections() == 1);
      CHECK(adj_list.connection_point_empty(0, 0) == false);
    }
    SUBCASE("test remove node") {
      int i = 0;
      AdjacencyList<int, FixedEdges> adj_list;
      adj_list.add_node(i, 1);
      adj_list.remove_node(0);
      CHECK(adj_list.get_num_nodes() == 0);
      CHECK_THROWS_AS(auto n = adj_list.get_node(0), GraphException);
      CHECK_THROWS_AS(auto d = adj_list.get_node_data(0), GraphException);
      CHECK_THROWS_AS(adj_list.remove_node(0), GraphException);
    }
    SUBCASE("test remove node 2") {
      int i = 0;
      AdjacencyList<int, FixedEdges> adj_list;
      adj_list.add_node(i, 1);
      adj_list.remove_node(0);
      i = 5;
      adj_list.add_node(i, 1);
      CHECK(adj_list.get_node_data(0) == 5);
      adj_list.add_node(i, 1);
      CHECK(adj_list.get_node_data(1) == 5);
      adj_list.remove_node(0);
      CHECK(adj_list.get_node_data(1) == 5);
      adj_list.add_node(i, 1);
      CHECK(adj_list.get_node_data(2) == 5);
    }
    SUBCASE("test remove connection") {
      int i = 0, j = 1;
      AdjacencyList<int, FixedEdges> adj_list;
      adj_list.add_node(i, 1);
      adj_list.add_node(j, 1);
      adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
      CHECK(adj_list.are_nodes_connected(0, 1));
      auto const &con1 = adj_list.get_node_connections(0);
      CHECK(con1.size() == 1);
      CHECK(con1[0] != nullptr);
      adj_list.remove_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
      CHECK(adj_list.are_nodes_connected(0, 1) == false);
      auto const &con2 = adj_list.get_node_connections(0);
      CHECK(con2.size() == 1);
      CHECK(con2[0] == nullptr);
    }
    SUBCASE("test remove node with a connection") {
      int i = 0, j = 1;
      AdjacencyList<int, FixedEdges> adj_list;
      adj_list.add_node(i, 1);
      adj_list.add_node(j, 1);
      adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
      adj_list.remove_node(0);
      CHECK_THROWS_AS(auto f = adj_list.connection_point_empty(0, 0),
                      GraphException);
      CHECK(adj_list.connection_point_empty(1, 0));
    }
  }
  SUBCASE("test dynamic edges") {
    int i = 0, j = 1;
    AdjacencyList<int, DynamicEdges> adj_list;
    adj_list.add_node(i, 0);
    adj_list.add_node(j, 0);
    adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
  }
  SUBCASE("check add_connection errors") {
    int i = 0, j = 1;
    AdjacencyList<int, FixedEdges> adj_list;
    adj_list.add_node(i, 2);
    adj_list.add_node(j, 2);
    adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
    // connection point 1 is filled
    CHECK_THROWS_AS(
        adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0}),
        GraphException);
    // connection point 2 is filled
    CHECK_THROWS_AS(
        adj_list.add_connection(ConnectionPoint{0, 1}, ConnectionPoint{1, 0}),
        GraphException);
    // node 2 does not exist
    CHECK_THROWS_AS(
        adj_list.add_connection(ConnectionPoint{2, 1}, ConnectionPoint{1, 1}),
        GraphException);
    // edge 2 does not exist
    CHECK_THROWS_AS(
        adj_list.add_connection(ConnectionPoint{1, 2}, ConnectionPoint{1, 1}),
        GraphException);
  }
  SUBCASE("test remove_connection errors") {
    int i = 0, j = 1;
    AdjacencyList<int, FixedEdges> adj_list;
    adj_list.add_node(i, 2);
    adj_list.add_node(j, 2);
    adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
    CHECK_THROWS_AS(adj_list.remove_connection(ConnectionPoint{0, 1},
                                               ConnectionPoint{1, 0}),
                    GraphException);
  }
  SUBCASE("test copying") {
    int i = 0, j = 1;
    FixedEdged_AL<int> adj_list;
    adj_list.add_node(i, 3);
    adj_list.add_node(j, 3);
    adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
    FixedEdged_AL<int> adj_list_2 = adj_list;
    adj_list.remove_node(0);
    CHECK(adj_list_2.get_num_connections() == 1);
    CHECK(adj_list_2.get_num_nodes() == 2);
  }

  /*


  SUBCASE("test directed") {
    auto adj_list = data_structure::FixedEdged_DAL<int>();
    adj_list.add_node(0, 3);
    adj_list.add_node(1, 3, 0, data_structure::NodeIndexandEdge{0, 0});

    CHECK(adj_list.get_num_nodes() == 2);

    CHECK(adj_list.has_parent(0) == false);
    CHECK(adj_list.has_parent(1) == true);
  }



  SUBCASE("test directed copy") {
    auto adj_list = data_structure::FixedEdged_DAL<int>();
    adj_list.add_node(0, 3);
    adj_list.add_node(1, 3, 0, data_structure::NodeIndexandEdge{0, 0});

    auto adj_list_2 = data_structure::FixedEdged_DAL<int>(adj_list);

    adj_list.remove_node(0);

    CHECK(adj_list_2.get_num_nodes() == 2);

    CHECK(adj_list_2.has_parent(0) == false);
    CHECK(adj_list_2.has_parent(1) == true);

    adj_list_2.add_node(2, 3, 0, data_structure::NodeIndexandEdge{1, 1});

    CHECK(adj_list_2.get_num_nodes() == 3);
  }  */

  SUBCASE("test directed copying") {
    int i = 0, j = 1;
    FixedEdged_DAL<int> adj_list;
    adj_list.add_node(i, 3);
    adj_list.add_node(j, 3, 0, ConnectionPoint{0, 0});
    FixedEdged_DAL<int> adj_list_2 = adj_list;
    adj_list.remove_node(0);
    CHECK(adj_list_2.get_num_connections() == 1);
    CHECK(adj_list_2.get_num_nodes() == 2);
  }
}