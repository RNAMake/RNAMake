//
// Created by Joseph Yesselman on 10/22/17.
//

#include "../common.hpp"
#include <iostream>

#include <data_structure/graph/graph_iter_list.h>

TEST_CASE("Test Graph Data Structure ") {
  using namespace data_structure::graph;
  SUBCASE("test trival case") {
    int i = 10, j = 10;
    FixedEdged_AL<int> adj_list;
    adj_list.add_node(i, 5);
    adj_list.add_node(j, 5);

    auto iter_list = IterList<int, FixedEdged_AL<int>>();
    iter_list.transversal(adj_list, 0);
    std::vector<Index> indexes;
    for (auto const &n : iter_list) {
      indexes.push_back(n->get_index());
    }
    CHECK(indexes[0] == 0);
    CHECK(indexes[1] == 1);
  }
  SUBCASE("test with undirected adjacency list") {
    std::vector<int> data = {0, 1, 2, 3, 4, 5};
    FixedEdged_AL<int> adj_list;
    adj_list.add_node(data[0], 5);
    adj_list.add_node(data[1], 5);
    adj_list.add_node(data[2], 5);
    adj_list.add_node(data[3], 5);
    adj_list.add_node(data[4], 5);
    adj_list.add_node(data[5], 5);

    adj_list.add_connection(ConnectionPoint{0, 0}, ConnectionPoint{1, 0});
    adj_list.add_connection(ConnectionPoint{1, 1}, ConnectionPoint{2, 0});
    adj_list.add_connection(ConnectionPoint{2, 1}, ConnectionPoint{3, 0});
    adj_list.add_connection(ConnectionPoint{3, 1}, ConnectionPoint{1, 2});
    adj_list.add_connection(ConnectionPoint{3, 2}, ConnectionPoint{4, 0});
    adj_list.add_connection(ConnectionPoint{4, 1}, ConnectionPoint{0, 2});

    /* adj_list connectivity
     *
     * 0 <-> 1 <-> 2 <-> 3 <-> 4   5
     * ^     ^           ^     ^
     * |     |           |     |
     * |     -------------     |
     * -------------------------
     */
    IterList<int, FixedEdged_AL<int>> iter_list;
    iter_list.transversal(adj_list, 0);
    auto count = 0;
    for (auto const &n : iter_list) {
      count += 1;
    }
    CHECK(count == adj_list.get_num_nodes());
    iter_list.path_transversal(adj_list, 0, 2);
    auto path = std::vector<Index>();
    auto target = std::vector<Index>{0, 1, 2};
    for (auto const &n : iter_list) {
      path.push_back(n->get_index());
    }
    CHECK(path == target);

    // no viable path
    REQUIRE_THROWS_AS(iter_list.path_transversal(adj_list, 0, 5),
                      GraphException);

    iter_list.path_transversal(adj_list, 0, 3);
    path = std::vector<Index>();
    target = std::vector<Index>{0, 1, 3};
    for (auto const &n : iter_list) {
      path.push_back(n->get_index());
    }
    CHECK(path == target);
  }
  SUBCASE("test with directed adjacency list") {
    std::vector<int> data = {0, 1, 2, 3, 4, 5};
    FixedEdged_DAL<int> adj_list;
    adj_list.add_node(data[0], 5);
    adj_list.add_node(data[1], 5, 0, ConnectionPoint{0, 0});
    adj_list.add_node(data[2], 5, 0, ConnectionPoint{0, 1});
    adj_list.add_node(data[3], 5, 0, ConnectionPoint{2, 1});
    adj_list.add_node(data[4], 5, 0, ConnectionPoint{3, 1});
    adj_list.add_node(data[5], 5, 0, ConnectionPoint{3, 2});
    adj_list.add_connection(ConnectionPoint{5, 1}, ConnectionPoint{1, 1});

    /* adj_list connectivity
     *
     *      5 -------
     *      ^       |
     *      |       |
     * 2 -> 3 -> 4  |
     * ^            |
     * |            |
     * 0 -> 1 -------
     *
     */
    DirectedIterList<int, FixedEdged_DAL<int>> iter_list;
    iter_list.transversal(adj_list, 0);
    auto path = std::vector<Index>();
    auto target = std::vector<Index>{0, 1, 2, 3, 4, 5};
    for (auto const &n : iter_list) {
      path.push_back(n->get_index());
    }
    CHECK(path == target);

    iter_list.path_transversal(adj_list, 0, 4);
    path = std::vector<Index>();
    target = std::vector<Index>{0, 2, 3, 4};
    for (auto const &n : iter_list) {
      path.push_back(n->get_index());
    }
    CHECK(path == target);

    iter_list.path_transversal(adj_list, 0, 5);
    path = std::vector<Index>();
    target = std::vector<Index>{0, 2, 3, 5};
    for (auto const &n : iter_list) {
      path.push_back(n->get_index());
    }
    CHECK(path == target);

    path = std::vector<Index>();
    target = std::vector<Index>{2, 3, 5, 4};
    iter_list.sub_graph_transversal(adj_list, 2, 5);
    for (auto const &n : iter_list) {
      path.push_back(n->get_index());
    }
    CHECK(path == target);
  }
}