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
    std::vector<int> data = {0, 1, 2, 3, 4, 5, 6, 7};
    auto g = FixedEdgeDirectedGraph<int>();
    g.add_node(data[0], 3);
    g.add_node(data[1], 3, 0, ConnectionPoint{0, 0});
    g.set_node_data(0, data[2]);
    CHECK(g.get_node_data(0) == data[2]);
    g.get_node_data(0) = data[3];
    CHECK(g.get_node_data(0) == data[3]);
    g[0] = data[4];
    CHECK(g[0] == data[4]);
    g.setup_transversal(0);
    for (auto const &n : g) {
      g[n->get_index()] = data[5];
      CHECK(n->get_data() == data[5]);
    }
  }
  SUBCASE("test new types") {
    struct X {
      int x;
      int y;
    };
    FixedEdgeDirectedGraph<X> g1;
    X x1 = {1, 1};
    g1.add_node(x1, 3);
    CHECK(g1.get_node_data(0).x == 1);
    g1.setup_transversal(0);
    int i = 0;
    for (auto const &n : g1) {
      i++;
    }
    CHECK(i == 1);
  }
  SUBCASE("test iteration forwarding") {
    struct GraphStruct {
    public:
      typedef
          typename FixedEdgeDirectedGraph<int>::const_iterator const_iterator;
      const_iterator begin() const {
        graph_.setup_transversal(0);
        return graph_.begin();
      }
      const_iterator end() const { return graph_.end(); }

    public:
      int add_node() {
        int i = 1;
        return graph_.add_node(i, 1);
      }

    private:
      mutable FixedEdgeDirectedGraph<int> graph_;
    };

    auto gs = GraphStruct();
    gs.add_node();
    gs.add_node();
    int i = 0;
    for (auto &n : gs) {
      i++;
    }
    CHECK(i == 2);
  }
}