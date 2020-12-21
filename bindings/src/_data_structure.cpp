// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// base includes
#include <base/application.hpp>
#include <base/backtrace.h>
#include <base/cl_option.h>
#include <base/command_line_parser.hpp>
#include <base/env_manager.h>
#include <base/exception.h>
#include <base/file_io.h>
#include <base/log.h>
#include <base/option.h>
#include <base/settings.h>
#include <base/string.h>
#include <base/sys_interface.h>
#include <base/vector_container.h>

// data_structure includes
#include <data_structure/graph.h>
#include <data_structure/graph_adjacency_list.h>
#include <data_structure/graph_base.h>
#include <data_structure/graph_iter_list.h>

// data_structure::graph includes
#include <data_structure/graph/graph.h>
#include <data_structure/graph/graph_node.fwd.h>
#include <data_structure/graph/graph_node.h>

// data_structure::tree includes
#include <data_structure/tree/tree.h>
#include <data_structure/tree/tree_node.h>

// eternabot includes
#include <eternabot/feature_generator.h>
#include <eternabot/scorer.h>
#include <eternabot/sequence_designer.h>
#include <eternabot/strategy.h>
#include <eternabot/strategy/a_basic_test.h>
#include <eternabot/strategy/berex_test.h>
#include <eternabot/strategy/clear_plot.h>
#include <eternabot/strategy/direction_of_gc.h>
#include <eternabot/strategy/num_of_yellow.h>

// io includes
#include <util/csv.h>

// io::detail includes
#include <util/csv.h>

// io::error includes
#include <util/csv.h>

// math includes
#include <math/euler.h>
#include <math/hashing.h>
#include <math/numerical.h>
#include <math/quaternion.h>
#include <math/stats.h>
#include <math/transform.h>
#include <math/xyz_matrix.h>
#include <math/xyz_vector.h>

using namespace data_structure;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_data_structure,m) {
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// data_structure
    ////////////////////////////////////////////////////////////////////////////////////////////////////
	// classes
/*
        py::class_<AdjacencyList, std::shared_ptr<AdjacencyList>>(m, "AdjacencyList")
		// ctors
		.def(py::init<>())
		.def(py::init<AdjacencyList const &>())
		// methods
		.def("begin",[] (AdjacencyList const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (AdjacencyList const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("add_node",[] (AdjacencyList  & ptr, DataType const & d, Size n_edges) -> Index {
		 return ptr.add_node(d, n_edges); } )
		.def("add_edge",[] (AdjacencyList  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
		ptr.add_edge(nie1, nie2); } )
		.def("remove_node",[] (AdjacencyList  & ptr, Index ni) {
		ptr.remove_node(ni); } )
		.def("remove_edge",[] (AdjacencyList  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
		ptr.remove_edge(nie1, nie2); } )
		.def("get_num_nodes",[] (AdjacencyList const & ptr) ->  size_t {
		 return ptr.get_num_nodes(); } )
		.def("get_num_edges",[] (AdjacencyList const & ptr) -> size_t {
		 return ptr.get_num_edges(); } )
		.def("get_node_edges",[] (AdjacencyList const & ptr, Index ni) -> std::vector<Edge const *> const & {
		 return ptr.get_node_edges(ni); } )
		.def("get_node_data",[] (AdjacencyList const & ptr, Index ni) -> DataType const & {
		 return ptr.get_node_data(ni); } )
		.def("get_node_data",[] (AdjacencyList  & ptr, Index ni) -> DataType & {
		 return ptr.get_node_data(ni); } )
		.def("get_node",[] (AdjacencyList const & ptr, Index ni) -> Node<DataType> const & {
		 return ptr.get_node(ni); } )
		.def("get_node",[] (AdjacencyList  & ptr, Index ni) -> Node<DataType> & {
		 return ptr.get_node(ni); } )
		.def("get_connected_node_info",[] (AdjacencyList const & ptr, NodeIndexandEdge const & nei) -> NodeIndexandEdge {
		 return ptr.get_connected_node_info(nei); } )
		.def("edge_between_nodes",[] (AdjacencyList const & ptr, Index n1, Index n2) -> bool {
		 return ptr.edge_between_nodes(n1, n2); } )
		.def("edge_index_empty",[] (AdjacencyList const & ptr, Index ni, Index ei) -> bool {
		 return ptr.edge_index_empty(ni, ei); } )
		// operators
		.def(py::self = py::self)
		;

        py::class_<DirectedAdjacencyList, std::shared_ptr<DirectedAdjacencyList>>(m, "DirectedAdjacencyList")
		// ctors
		.def(py::init<>())
		.def(py::init<DirectedAdjacencyList const &>())
		// methods
		.def("add_node",[] (DirectedAdjacencyList  & ptr, DataType const & d, Size n_edges) -> Index {
		 return ptr.add_node(d, n_edges); } )
		.def("add_node",[] (DirectedAdjacencyList  & ptr, DataType const & d, Size n_edges, Index n_end_index, NodeIndexandEdge const & pie) -> Index {
		 return ptr.add_node(d, n_edges, n_end_index, pie); } )
		.def("remove_node",[] (DirectedAdjacencyList  & ptr, Index ni) {
		ptr.remove_node(ni); } )
		.def("has_parent",[] (DirectedAdjacencyList const & ptr, Index n_index) -> bool {
		 return ptr.has_parent(n_index); } )
		.def("get_parent_index",[] (DirectedAdjacencyList const & ptr, Index n_index) -> Index {
		 return ptr.get_parent_index(n_index); } )
		.def("get_parent_end_index",[] (DirectedAdjacencyList const & ptr, Index n_index) -> Index {
		 return ptr.get_parent_end_index(n_index); } )
		// operators
		.def(py::self = py::self)
		;

        py::class_<DirectedGraph, std::shared_ptr<DirectedGraph>>(m, "DirectedGraph")
		// ctors
		.def(py::init<>())
		.def(py::init<DirectedGraph const &>())
		// methods
		.def("setup_sub_graph_transversal",[] (DirectedGraph  & ptr, Index start_n, Index end_n) {
		ptr.setup_sub_graph_transversal(start_n, end_n); } )
		.def("add_node",[] (DirectedGraph  & ptr, DataType const & d, Size n_edges, Index n_end_index, NodeIndexandEdge const & pie) ->  Index {
		 return ptr.add_node(d, n_edges, n_end_index, pie); } )
		.def("add_node",[] (DirectedGraph  & ptr, DataType const & d, Size n_edges) ->  Index {
		 return ptr.add_node(d, n_edges); } )
		.def("has_parent",[] (DirectedGraph const & ptr, Index ni) ->  bool {
		 return ptr.has_parent(ni); } )
		.def("get_parent_index",[] (DirectedGraph const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_index(ni); } )
		.def("get_parent_end_index",[] (DirectedGraph const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_end_index(ni); } )
		.def("get_root_indexes",[] (DirectedGraph  & ptr) -> Indexes {
		 return ptr.get_root_indexes(); } )
		;

        py::class_<DirectedIterList, std::shared_ptr<DirectedIterList>>(m, "DirectedIterList")
		// ctors
		.def(py::init<>())
		// methods
		.def("transversal",[] (DirectedIterList  & ptr, AdjacencyListType & adj_list, Index start_n) {
		ptr.transversal(adj_list, start_n); } )
		.def("sub_graph_transversal",[] (DirectedIterList  & ptr, AdjacencyListType & adj_list, Index start_n, Index end_n) {
		ptr.sub_graph_transversal(adj_list, start_n, end_n); } )
		;
*/
        py::class_<DynamicEdges, std::shared_ptr<DynamicEdges>>(m, "DynamicEdges")
		;

        py::class_<Edge, std::shared_ptr<Edge>>(m, "Edge")
		// ctors
		.def(py::init<Index,Index,Index,Index>())
		// methods
		.def("partner",[] (Edge const & ptr, Index index) -> Index {
		 return ptr.partner(index); } )
		.def("end_index",[] (Edge const & ptr, Index index) -> Index {
		 return ptr.end_index(index); } )
		.def("to_str",[] (Edge const & ptr) -> String {
		 return ptr.to_str(); } )
		// operators
		.def(py::self == py::self)
		// public attributes
		.def_readwrite("node_i", &Edge::node_i)
		.def_readwrite("node_j", &Edge::node_j)
		.def_readwrite("edge_i", &Edge::edge_i)
		.def_readwrite("edge_j", &Edge::edge_j)
		;

        py::class_<FixedEdges, std::shared_ptr<FixedEdges>>(m, "FixedEdges")
		;
/*
        py::class_<IterList, std::shared_ptr<IterList>>(m, "IterList")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (IterList  & ptr) -> iterator {
		 return ptr.begin(); } )
		.def("end",[] (IterList  & ptr) -> iterator {
		 return ptr.end(); } )
		.def("begin",[] (IterList const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (IterList const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("transversal",[] (IterList  & ptr, AdjacencyListType & adj_list, Index start_n) {
		ptr.transversal(adj_list, start_n); } )
		.def("path_transversal",[] (IterList  & ptr, AdjacencyListType & adj_list, Index start_n, Index end_n) {
		ptr.path_transversal(adj_list, start_n, end_n); } )
		;

        py::class_<Node, std::shared_ptr<Node>>(m, "Node")
		// ctors
		.def(py::init<DataType const &,Index const>())
		// methods
		.def("data",[] (Node const & ptr) ->  DataType const & {
		 return ptr.data(); } )
		.def("data",[] (Node  & ptr) ->  DataType & {
		 return ptr.data(); } )
		.def("index",[] (Node const & ptr) ->  Index {
		 return ptr.index(); } )
		;
*/
        py::class_<NodeIndexandEdge, std::shared_ptr<NodeIndexandEdge>>(m, "NodeIndexandEdge")
		// methods
		.def("to_str",[] (NodeIndexandEdge const & ptr) -> String {
		 return ptr.to_str(); } )
		// operators
		.def(py::self == py::self)
		// public attributes
		.def_readwrite("node_index", &NodeIndexandEdge::node_index)
		.def_readwrite("edge_index", &NodeIndexandEdge::edge_index)
		;
/*
        py::class_<NodeIndexandEdgeCompare, std::shared_ptr<NodeIndexandEdgeCompare>>(m, "NodeIndexandEdgeCompare")
		// operators
		.def(py::self () data_structure::NodeIndexandEdge const &)
		;

        py::class_<UndirectedGraph, std::shared_ptr<UndirectedGraph>>(m, "UndirectedGraph")
		// ctors
		.def(py::init<>())
		.def(py::init<UndirectedGraph const &>())
		;

        py::class_<VisitedNode, std::shared_ptr<VisitedNode>>(m, "VisitedNode")
		// methods
		.def("index_in_path",[] (VisitedNode  & ptr, Index i) -> bool {
		 return ptr.index_in_path(i); } )
		.def("path_length",[] (VisitedNode  & ptr) -> int {
		 return ptr.path_length(); } )
		// public attributes
		.def_readwrite("parent", &VisitedNode::parent)
		.def_readwrite("index", &VisitedNode::index)
		;

        py::class_<_Graph, std::shared_ptr<_Graph>>(m, "_Graph")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (_Graph  & ptr) -> iterator {
		 return ptr.begin(); } )
		.def("end",[] (_Graph  & ptr) -> iterator {
		 return ptr.end(); } )
		.def("begin",[] (_Graph const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (_Graph const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("setup_transversal",[] (_Graph  & ptr, Index start_n) {
		ptr.setup_transversal(start_n); } )
		.def("setup_path_transversal",[] (_Graph  & ptr, Index start_n, Index end_n) {
		ptr.setup_path_transversal(start_n, end_n); } )
		.def("add_node",[] (_Graph  & ptr, DataType const & d, Size n_edges) ->  Index {
		 return ptr.add_node(d, n_edges); } )
		.def("add_edge",[] (_Graph  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
		ptr.add_edge(nie1, nie2); } )
		.def("remove_node",[] (_Graph  & ptr, Index ni) {
		ptr.remove_node(ni); } )
		.def("remove_edge",[] (_Graph  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
		ptr.remove_edge(nie1, nie2); } )
		.def("get_num_nodes",[] (_Graph const & ptr) ->  size_t {
		 return ptr.get_num_nodes(); } )
		.def("get_num_edges",[] (_Graph const & ptr) ->  size_t {
		 return ptr.get_num_edges(); } )
		.def("get_node_edges",[] (_Graph const & ptr, Index ni) ->  std::vector<Edge const *> const & {
		 return ptr.get_node_edges(ni); } )
		.def("get_node",[] (_Graph const & ptr, Index ni) ->  Node<DataType> const & {
		 return ptr.get_node(ni); } )
		.def("get_node_data",[] (_Graph const & ptr, Index ni) ->  DataType const & {
		 return ptr.get_node_data(ni); } )
		.def("get_node_data",[] (_Graph  & ptr, Index ni) ->  DataType & {
		 return ptr.get_node_data(ni); } )
		.def("get_connected_node_info",[] (_Graph const & ptr, NodeIndexandEdge const & nei) ->  NodeIndexandEdge {
		 return ptr.get_connected_node_info(nei); } )
		.def("edge_between_nodes",[] (_Graph const & ptr, Index n1, Index n2) ->  bool {
		 return ptr.edge_between_nodes(n1, n2); } )
		.def("edge_index_empty",[] (_Graph const & ptr, Index ni, Index ei) ->  bool {
		 return ptr.edge_index_empty(ni, ei); } )
		;
*/

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// data_structure::graph
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes
/*
        py::class_<Graph, std::shared_ptr<Graph>>(m, "Graph")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (Graph  & ptr) -> iterator {
		 return ptr.begin(); } )
		.def("end",[] (Graph  & ptr) -> iterator {
		 return ptr.end(); } )
		.def("begin",[] (Graph const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (Graph const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("transverse",[] (Graph const & ptr, data_structure::graphGraphNodeOP<DataType> const &) -> iterator {
		 return ptr.transverse(&); } )
		.def("size",[] (Graph const & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("get_node",[] (Graph const & ptr, int index) ->  GraphNodeOP<DataType> const & {
		 return ptr.get_node(index); } )
		.def("oldest_node",[] (Graph  & ptr) ->  GraphNodeOP<DataType> {
		 return ptr.oldest_node(); } )
		.def("increase_level",[] (Graph  & ptr) {
		ptr.increase_level(); } )
		.def("decrease_level",[] (Graph  & ptr) {
		ptr.decrease_level(); } )
		.def("nodes",[] (Graph const & ptr) ->  GraphNodeOPs<DataType> const & {
		 return ptr.nodes(); } )
		.def("connections",[] (Graph const & ptr) ->  GraphConnectionOPs<DataType> const & {
		 return ptr.connections(); } )
		.def("last_node",[] (Graph  & ptr) ->  GraphNodeOP<DataType> const & {
		 return ptr.last_node(); } )
		.def("level",[] (Graph  & ptr) ->  int {
		 return ptr.level(); } )
		.def("index",[] (Graph  & ptr) ->  int {
		 return ptr.index(); } )
		.def("index",[] (Graph  & ptr, int nindex) {
		ptr.index(nindex); } )
		;

        py::class_<GraphConnection, std::shared_ptr<GraphConnection>>(m, "GraphConnection")
		// ctors
		.def(py::init<GraphNodeOP<DataType>,GraphNodeOP<DataType>,int,int>())
		// methods
		.def("disconnect",[] (GraphConnection  & ptr) {
		ptr.disconnect(); } )
		.def("partner",[] (GraphConnection  & ptr, int i) -> GraphNodeOP<DataType> const & {
		 return ptr.partner(i); } )
		.def("end_index",[] (GraphConnection  & ptr, int n_index) -> int {
		 return ptr.end_index(n_index); } )
		.def("node_1",[] (GraphConnection  & ptr) ->  GraphNodeOP<DataType> const & {
		 return ptr.node_1(); } )
		.def("node_2",[] (GraphConnection  & ptr) ->  GraphNodeOP<DataType> const & {
		 return ptr.node_2(); } )
		.def("end_index_1",[] (GraphConnection  & ptr) ->  int {
		 return ptr.end_index_1(); } )
		.def("end_index_2",[] (GraphConnection  & ptr) ->  int {
		 return ptr.end_index_2(); } )
		;

        py::class_<GraphDynamic, std::shared_ptr<GraphDynamic>>(m, "GraphDynamic")
		// ctors
		.def(py::init<>())
		// methods
		.def("add_data",[] (GraphDynamic  & ptr, DataType const & data, int parent_index, int orphan) ->  int {
		 return ptr.add_data(data, parent_index, orphan); } )
		.def("connect",[] (GraphDynamic  & ptr, int i, int j) {
		ptr.connect(i, j); } )
		;

        py::class_<GraphIterator, std::shared_ptr<GraphIterator>>(m, "GraphIterator")
		// ctors
		.def(py::init<>())
		// operators
		.def(py::self ++ int)
		.def(py::self == py::self)
		.def(py::self != py::self)
		;

        py::class_<GraphNode, std::shared_ptr<GraphNode>>(m, "GraphNode")
		// ctors
		.def(py::init<int,int,size_t>())
		.def(py::init<DataType const &,int,int,size_t>())
		// methods
		.def("add_connection",[] (GraphNode  & ptr, data_structure::graphGraphConnectionOP<DataType> const &, int pos) {
		ptr.add_connection(&, pos); } )
		.def("remove_connection",[] (GraphNode  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection) {
		ptr.remove_connection(connection); } )
		.def("available_children_pos",[] (GraphNode const & ptr) ->  Ints {
		 return ptr.available_children_pos(); } )
		.def("available_pos",[] (GraphNode  & ptr, int pos) ->  int {
		 return ptr.available_pos(pos); } )
		.def("parent",[] (GraphNode  & ptr) ->  GraphNodeOP<DataType> {
		 return ptr.parent(); } )
		.def("parent_index",[] (GraphNode  & ptr) ->  int {
		 return ptr.parent_index(); } )
		.def("parent_end_index",[] (GraphNode  & ptr) ->  int {
		 return ptr.parent_end_index(); } )
		.def("unset_connections",[] (GraphNode  & ptr) {
		ptr.unset_connections(); } )
		.def("connected",[] (GraphNode  & ptr, data_structure::graphGraphNodeOP<DataType> const & n) ->  GraphConnectionOP<DataType> {
		 return ptr.connected(n); } )
		.def("index",[] (GraphNode const & ptr) ->  int {
		 return ptr.index(); } )
		.def("level",[] (GraphNode const & ptr) ->  int {
		 return ptr.level(); } )
		.def("data",[] (GraphNode const & ptr) ->  DataType const & {
		 return ptr.data(); } )
		.def("data",[] (GraphNode  & ptr) ->  DataType & {
		 return ptr.data(); } )
		.def("connections",[] (GraphNode const & ptr) ->  GraphConnectionOPs<DataType> const & {
		 return ptr.connections(); } )
		.def("index",[] (GraphNode  & ptr, int index) {
		ptr.index(index); } )
		;

        py::class_<GraphNodeCompare, std::shared_ptr<GraphNodeCompare>>(m, "GraphNodeCompare")
		// operators
		.def(py::self () GraphNodeOP<DataType> const &)
		;

        py::class_<GraphNodeDynamic, std::shared_ptr<GraphNodeDynamic>>(m, "GraphNodeDynamic")
		// ctors
		.def(py::init<DataType const &,int,int>())
		// methods
		.def("add_connection",[] (GraphNodeDynamic  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection, int pos) {
		ptr.add_connection(connection, pos); } )
		.def("remove_connection",[] (GraphNodeDynamic  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection) {
		ptr.remove_connection(connection); } )
		;

        py::class_<GraphNodeStatic, std::shared_ptr<GraphNodeStatic>>(m, "GraphNodeStatic")
		// ctors
		.def(py::init<DataType const &,int,int,int>())
		.def(py::init<GraphNode<DataType> const &>())
		// methods
		.def("add_connection",[] (GraphNodeStatic  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection, int pos) {
		ptr.add_connection(connection, pos); } )
		.def("remove_connection",[] (GraphNodeStatic  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection) {
		ptr.remove_connection(connection); } )
		;

        py::class_<GraphNodeType, std::shared_ptr<GraphNodeType>>(m, "GraphNodeType")
		;

        py::class_<GraphStatic, std::shared_ptr<GraphStatic>>(m, "GraphStatic")
		// ctors
		.def(py::init<>())
		.def(py::init<GraphStatic<DataType> const &>())
		// methods
		.def("add_data",[] (GraphStatic  & ptr, DataType const & data, int parent_index, int parent_pos, int child_pos, int n_children, int orphan, int index) ->  int {
		 return ptr.add_data(data, parent_index, parent_pos, child_pos, n_children, orphan, index); } )
		.def("connect",[] (GraphStatic  & ptr, int i, int j, int i_pos, int j_pos) {
		ptr.connect(i, j, i_pos, j_pos); } )
		.def("check_pos_is_valid",[] (GraphStatic  & ptr, data_structure::graphGraphNodeOP<DataType> const & n, int & pos) ->  int {
		 return ptr.check_pos_is_valid(n, pos); } )
		.def("get_available_pos",[] (GraphStatic  & ptr, data_structure::graphGraphNodeOP<DataType> const & n, int & pos) ->  Ints {
		 return ptr.get_available_pos(n, pos); } )
		.def("remove_node",[] (GraphStatic  & ptr, int pos) {
		ptr.remove_node(pos); } )
		.def("remove_level",[] (GraphStatic  & ptr, int level) {
		ptr.remove_level(level); } )
		;


*/

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// data_structure::tree
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes
/*
        py::class_<Tree, std::shared_ptr<Tree>>(m, "Tree")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (Tree  & ptr) -> iterator {
		 return ptr.begin(); } )
		.def("end",[] (Tree  & ptr) -> iterator {
		 return ptr.end(); } )
		.def("begin",[] (Tree const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (Tree const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("get_node",[] (Tree  & ptr, int index) ->  TreeNodeOP<DataType> const & {
		 return ptr.get_node(index); } )
		.def("remove_node",[] (Tree  & ptr, data_structure::treeTreeNodeOP<DataType> const & n) {
		ptr.remove_node(n); } )
		.def("remove_node",[] (Tree  & ptr, int pos) {
		ptr.remove_node(pos); } )
		.def("remove_level",[] (Tree  & ptr, int level) {
		ptr.remove_level(level); } )
		.def("increase_level",[] (Tree  & ptr) {
		ptr.increase_level(); } )
		.def("decrease_level",[] (Tree  & ptr) {
		ptr.decrease_level(); } )
		.def("size",[] (Tree  & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("last_node",[] (Tree  & ptr) ->  TreeNodeOP<DataType> const & {
		 return ptr.last_node(); } )
		.def("level",[] (Tree  & ptr) ->  int {
		 return ptr.level(); } )
		;

        py::class_<TreeDynamic, std::shared_ptr<TreeDynamic>>(m, "TreeDynamic")
		// ctors
		.def(py::init<>())
		// methods
		.def("add_data",[] (TreeDynamic  & ptr, DataType const & data, int parent_index) ->  int {
		 return ptr.add_data(data, parent_index); } )
		;

        py::class_<TreeNode, std::shared_ptr<TreeNode>>(m, "TreeNode")
		// ctors
		.def(py::init<DataType const &,int,int,size_t>())
		// methods
		.def("add_child",[] (TreeNode  & ptr, data_structure::treeTreeNodeOP<DataType> const &, int pos) {
		ptr.add_child(&, pos); } )
		.def("remove_child",[] (TreeNode  & ptr, data_structure::treeTreeNodeOP<DataType> const &) {
		ptr.remove_child(&); } )
		.def("leaf",[] (TreeNode  & ptr) -> bool {
		 return ptr.leaf(); } )
		.def("available_children_pos",[] (TreeNode const & ptr) ->  Ints {
		 return ptr.available_children_pos(); } )
		.def("available_pos",[] (TreeNode  & ptr, int pos) ->  int {
		 return ptr.available_pos(pos); } )
		.def("unset_children",[] (TreeNode  & ptr) {
		ptr.unset_children(); } )
		.def("index",[] (TreeNode const & ptr) ->  int {
		 return ptr.index(); } )
		.def("level",[] (TreeNode const & ptr) ->  int {
		 return ptr.level(); } )
		.def("data",[] (TreeNode const & ptr) ->  DataType const & {
		 return ptr.data(); } )
		.def("data",[] (TreeNode  & ptr) ->  DataType & {
		 return ptr.data(); } )
		.def("children",[] (TreeNode  & ptr) ->  TreeNodeOPs<DataType> const & {
		 return ptr.children(); } )
		.def("parent",[] (TreeNode  & ptr) ->  TreeNodeOP<DataType> const & {
		 return ptr.parent(); } )
		.def("parent_index",[] (TreeNode  & ptr) ->  int {
		 return ptr.parent_index(); } )
		.def("parent_end_index",[] (TreeNode  & ptr) ->  int {
		 return ptr.parent_end_index(); } )
		.def("parent",[] (TreeNode  & ptr, data_structure::treeTreeNodeOP<DataType> const & p) {
		ptr.parent(p); } )
		.def("index",[] (TreeNode  & ptr, int index) {
		ptr.index(index); } )
		;

        py::class_<TreeNodeDynamic, std::shared_ptr<TreeNodeDynamic>>(m, "TreeNodeDynamic")
		// ctors
		.def(py::init<DataType const &,int,int>())
		// methods
		.def("add_child",[] (TreeNodeDynamic  & ptr, data_structure::treeTreeNodeOP<DataType> const & c, int pos) {
		ptr.add_child(c, pos); } )
		.def("remove_child",[] (TreeNodeDynamic  & ptr, data_structure::treeTreeNodeOP<DataType> const & child) {
		ptr.remove_child(child); } )
		.def("leaf",[] (TreeNodeDynamic  & ptr) ->  bool {
		 return ptr.leaf(); } )
		;

        py::class_<TreeNodeStatic, std::shared_ptr<TreeNodeStatic>>(m, "TreeNodeStatic")
		// ctors
		.def(py::init<DataType const &,int,int,int>())
		// methods
		.def("add_child",[] (TreeNodeStatic  & ptr, data_structure::treeTreeNodeOP<DataType> const & child, int pos) {
		ptr.add_child(child, pos); } )
		.def("remove_child",[] (TreeNodeStatic  & ptr, data_structure::treeTreeNodeOP<DataType> const & child) {
		ptr.remove_child(child); } )
		.def("leaf",[] (TreeNodeStatic  & ptr) ->  bool {
		 return ptr.leaf(); } )
		;

        py::class_<TreeStatic, std::shared_ptr<TreeStatic>>(m, "TreeStatic")
		// ctors
		.def(py::init<>())
		.def(py::init<TreeStatic<DataType> const &>())
		// methods
		.def("add_data",[] (TreeStatic  & ptr, DataType const & data, int n_children, int parent_index, int parent_child_index) ->  int {
		 return ptr.add_data(data, n_children, parent_index, parent_child_index); } )
		.def("get_available_pos",[] (TreeStatic  & ptr, data_structure::treeTreeNodeOP<DataType> const & n, int pos) ->  Ints {
		 return ptr.get_available_pos(n, pos); } )
		.def("index",[] (TreeStatic  & ptr, int n_index) {
		ptr.index(n_index); } )
		;

*/

}