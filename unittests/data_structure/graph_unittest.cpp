//
// Created by Joseph Yesselman on 10/22/17.
//

#include "../common.hpp"

#include <data_structure/graph/graph.h>

#include <iostream>

TEST_CASE("Test Graph Data Structure ") {
  using namespace data_structure::graph;
  SUBCASE("test basic functionality") {
    int i = 0, j = 1;
    auto g = FixedEdgeDirectedGraph<int>();
    g.add_node(i, 3);
    g.add_node(j, 3, 0, ConnectionPoint{0, 0});

    CHECK(g.get_num_nodes() == 2);
    CHECK(g.has_parent(0) == false);
    CHECK(g.has_parent(1) == true);

    g.setup_transversal(0);
    auto path = std::vector<Index>();
    for (auto const &n : g) {
      path.push_back(n->get_index());
    }
    CHECK(path.size() == 2);

    auto g2 = FixedEdgeDirectedGraph<int>(g);
    g.remove_node(0);

    CHECK(g2.get_num_nodes() == 2);
    CHECK(g2.get_num_connections() == 1);

    g2.remove_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});

    CHECK(g2.get_num_connections() == 0);

    // test undirected graph
    i = 0, j = 1;
    auto g3 = FixedEdgeUndirectedGraph<int>();
    g3.add_node(i, 3);
    g3.add_node(j, 3);
    g3.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
  }
  SUBCASE("test setting new values") {
    int i = 0, j = 1, k = 2, l = 3;
    auto g = FixedEdgeDirectedGraph<int>();
    g.add_node(i, 3);
    g.add_node(j, 3, 0, ConnectionPoint{0, 0});
    g.set_node_data(0, k);
    CHECK(g.get_node_data(0) == k);
    g.get_node_data(0) = l;
    CHECK(g.get_node_data(0) == l);

  }

  /*SUBCASE("test new types") {

    struct X {
      int x;
      int y;
    };

    auto g1 = data_structure::FixedEdgeDirectedGraph<X>();
    g1.add_node(X{0, 0}, 3);

    CHECK(g1.get_node_data(0).x == 0);

    g1.setup_transversal(0);
    int i = 0;
    for (auto const &n : g1) {
      i++;
    }
  }

  SUBCASE("test iteration forwarding") {
    struct GraphStruct {
    public:
      GraphStruct() : graph_(data_structure::FixedEdgeDirectedGraph<int>()) {}

    public:
      typedef
          typename data_structure::FixedEdgeDirectedGraph<int>::const_iterator
              const_iterator;
      typedef typename data_structure::FixedEdgeDirectedGraph<int>::iterator
          iterator;

      iterator begin() { return graph_.begin(); }
      iterator end() { return graph_.end(); }

      const_iterator begin() const { return graph_.begin(); }
      const_iterator end() const { return graph_.end(); }

    private:
      data_structure::FixedEdgeDirectedGraph<int> graph_;
    };

    int i = 0;
    auto gs = GraphStruct();
    for (auto &n : gs) {
      i++;
    }
  }*/
}