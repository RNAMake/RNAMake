//
// Created by Joe Yesselman on 6/18/22.
//

#include "../common.hpp"
#include <iostream>

#include <data_structure/graph/graph_base.h>

TEST_CASE("test graph components") {
  using namespace data_structure::graph;
  SUBCASE("test node") {
    int num = 10;
    int num2 = 11;
    Node<int> n(num, 0);
    CHECK(n.get_data() == num);
    // can use non const setting
    n.set_data(num2);
    CHECK(n.get_data() != num);
  }
  SUBCASE("test connection") {
    SUBCASE("test proper initiation") {
      Connection c1 = {1, 2, 3, 4};
      CHECK(c1.node_i == 1);
      CHECK(c1.node_j == 2);
      CHECK(c1.edge_i == 3);
      CHECK(c1.edge_j == 4);
    }
    SUBCASE("check == ") {
      Connection c1 = {1, 2, 3, 4};
      Connection c2 = {1, 2, 3, 4};
      CHECK(c1 == c2);
      Connection c3 = {2, 2, 3, 4};
      CHECK(!(c1 == c3));
    }
    SUBCASE("test partner") {
      Connection c1 = {1, 2, 3, 4};
      CHECK(c1.get_partner(1) == 2);
      CHECK(c1.get_partner(2) == 1);
      CHECK_THROWS_AS(Index i = c1.get_partner(-1), GraphException);
    }
    SUBCASE("test get edge index") {
      Connection c1 = {1, 2, 3, 4};
      CHECK(c1.get_edge_index(1) == 3);
      CHECK(c1.get_edge_index(2) == 4);
      CHECK_THROWS_AS(Index i = c1.get_edge_index(-1), GraphException);
    }
  }

}