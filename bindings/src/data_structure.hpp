#ifndef PYBIND11_DATA_STRUCTURE_HPP
#define PYBIND11_DATA_STRUCTURE_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

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

String
capitalize_first(String const & input) {
    auto result = input;
    if (!std::isupper(result[0]))  {
        result[0] = std::toupper(result[0]);
    }
    return result;
}

namespace data_structure {
    namespace  py = pybind11 ;
    void
    add_bindings(py::module_&);

    template<
            typename DataType,
            typename AdjList = AdjacencyList<DataType, FixedEdges>,
            typename DirAdjList = DirectedAdjacencyList<DataType,FixedEdges>,
            typename DirGraph =  DirectedGraph<DataType,FixedEdges>,
            typename DirIterList = DirectedIterList<DataType, DirAdjList>,
            typename IteratorList = IterList<DataType, data_structure::FixedEdged_AL<DataType>>,
            typename NodeType = Node<DataType>,
            typename UnDirGraph = UndirectedGraph<DataType,FixedEdges>,
            typename Visited = typename IterList<DataType, data_structure::FixedEdged_AL<DataType>>::VisitedNode,
            typename Graph = _Graph<DataType, DirAdjList, IteratorList>
            >
    void
    add_graph_specialization(py::module_ & m , String const&  PyType) {
        //using AdjList = AdjacencyList<DataType, FixedEdges>;
        //using DirAdjList = DirectedAdjacencyList<DataType, FixedEdges>;
        //using DirGraph = DirectedGraph<int,FixedEdges>

        py::class_<NodeType, std::shared_ptr<NodeType>>(m, String(String("Node") + PyType).c_str())
                // ctors
                .def(py::init<DataType const &, Index const>(), py::arg("data"), py::arg("index"))
                        // methods
                .def("data",[] (NodeType const & ptr) ->  DataType const & { return ptr.data(); } )
                .def("data",[] (NodeType & ptr) ->  DataType & { return ptr.data(); } )
                .def("index",[] (NodeType const & ptr) ->  Index { return ptr.index(); } )
                ;

        py::class_<IteratorList, std::shared_ptr<IteratorList>>(m, String(String("IteratorList") + PyType).c_str())
                // ctors
                .def(py::init<>())
                        // methods
                .def("begin",[] (IteratorList  & ptr) {
                    return ptr.begin(); } )
                .def("end",[] (IteratorList  & ptr) {
                    return ptr.end(); } )
                .def("begin",[] (IteratorList const & ptr) {
                    return ptr.begin(); } )
                .def("end",[] (IteratorList const & ptr) {
                    return ptr.end(); } )
                .def("transversal",[] (IteratorList & ptr, AdjList & adj_list, Index start_n) {
                    ptr.transversal(adj_list, start_n); }, py::arg("adj_list"), py::arg("start_n") )
                .def("path_transversal",[] (IteratorList  & ptr, AdjList & adj_list, Index start_n, Index end_n) {
                    ptr.path_transversal(adj_list, start_n, end_n); }, py::arg("adj_list"), py::arg("start_n"), py::arg("end_n") )
                ;
        py::class_<DirIterList>(m, String(String("DirectedIterList") + PyType).c_str())
                // ctors
                .def(py::init<>())
                        // methods
                .def("transversal",[] (DirIterList  & ptr,  DirAdjList & adj_list, Index start_n) {
                    ptr.transversal(adj_list, start_n); }, py::arg("adj_list"), py::arg("start_n"))
                .def("sub_graph_transversal",[] (DirIterList  & ptr, DirAdjList & adj_list, Index start_n, Index end_n) {
                    ptr.sub_graph_transversal(adj_list, start_n, end_n); }, py::arg("adj_list"), py::arg("start_n"), py::arg("end_n") )
                ;

        py::class_<DirGraph, std::shared_ptr<DirGraph>>(m, String(String("DirectedGraph") + PyType).c_str())
                // ctors
                .def(py::init<>())
                .def(py::init<DirGraph const &>(), py::arg("g"))
                        // methods
                .def("setup_sub_graph_transversal",[] (DirGraph  & ptr, Index start_n, Index end_n) {
                    ptr.setup_sub_graph_transversal(start_n, end_n); }, py::arg("start_n"), py::arg("end_n") )
                .def("add_node",[] (DirGraph  & ptr, DataType const & d, Size n_edges, Index n_end_index, NodeIndexandEdge const & pie) ->  Index {
                    return ptr.add_node(d, n_edges, n_end_index, pie); }, py::arg("d"), py::arg("n_edges"), py::arg("n_end_index"), py::arg("pie") )
                .def("add_node",[] (DirGraph  & ptr, DataType const & d, Size n_edges) ->  Index {
                    return ptr.add_node(d, n_edges); }, py::arg("d"), py::arg("n_edges") )
                .def("has_parent",[] (DirGraph const & ptr, Index ni) ->  bool {
                    return ptr.has_parent(ni); }, py::arg("ni") )
                .def("get_parent_index",[] (DirGraph const & ptr, Index ni) ->  Index {
                    return ptr.get_parent_index(ni); }, py::arg("ni") )
                .def("get_parent_end_index",[] (DirGraph const & ptr, Index ni) ->  Index {
                    return ptr.get_parent_end_index(ni); }, py::arg("ni") )
                .def("get_root_indexes",[] (DirGraph  & ptr) -> Indexes {
                    return ptr.get_root_indexes(); } )
                ;

        py::class_<DirAdjList, std::shared_ptr<DirAdjList>>(m, String(String("DirectedAdjacencyList") +  PyType).c_str())
                // ctors
                .def(py::init<>())
                .def(py::init<DirAdjList const &>(), py::arg("abj_list"))
                        // methods
                .def("add_node",[] (DirAdjList  & ptr, DataType const & d, Size n_edges) -> Index {
                    return ptr.add_node(d, n_edges); }, py::arg("d"), py::arg("n_edges") )
                .def("add_node",[] (DirAdjList  & ptr, DataType const & d, Size n_edges, Index n_end_index, NodeIndexandEdge const & pie) -> Index {
                    return ptr.add_node(d, n_edges, n_end_index, pie); }, py::arg("d"), py::arg("n_edges"), py::arg("n_end_index"), py::arg("pie") )
                .def("remove_node",[] (DirAdjList  & ptr, Index ni) {
                    ptr.remove_node(ni); }, py::arg("ni") )
                .def("has_parent",[] (DirAdjList const & ptr, Index n_index) -> bool {
                    return ptr.has_parent(n_index); }, py::arg("n_index") )
                .def("get_parent_index",[] (DirAdjList const & ptr, Index n_index) -> Index {
                    return ptr.get_parent_index(n_index); }, py::arg("n_index") )
                .def("get_parent_end_index",[] (DirAdjList const & ptr, Index n_index) -> Index {
                    return ptr.get_parent_end_index(n_index); }, py::arg("n_index") )
            // operators
            // .def(py::self = py::self)
                ;

        py::class_<AdjList>(m, String(String("AdjacencyList") + PyType).c_str())
                // ctors
                .def(py::init<>())
                .def(py::init<AdjList const &>(), py::arg("adj_list"))
                        // methods
                .def("begin",[] (AdjList const & ptr)  { return ptr.begin(); } )
                .def("end",[] (AdjList const & ptr)  { return ptr.end(); } )
                .def("add_node",[] (AdjList  & ptr, DataType const & d, Size n_edges) -> Index { return ptr.add_node(d, n_edges); },
                     py::arg("d"), py::arg("n_edges")
                )
                .def("add_edge",[] (AdjList  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) { ptr.add_edge(nie1, nie2); },
                     py::arg("nie1"), py::arg("nie2")
                )
                .def("remove_node",[] (AdjList  & ptr, Index ni) {
                    ptr.remove_node(ni); }, py::arg("ni") )
                .def("remove_edge",[] (AdjList  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
                    ptr.remove_edge(nie1, nie2); }, py::arg("nie1"), py::arg("nie2") )
                .def("get_num_nodes",[] (AdjList const & ptr) ->  size_t { return ptr.get_num_nodes(); } )
                .def("get_num_edges",[] (AdjList const & ptr) -> size_t { return ptr.get_num_edges(); } )
                .def("get_node_edges",[] (AdjList const & ptr, Index ni) -> std::vector<Edge const *> const & {
                    return ptr.get_node_edges(ni); }, py::arg("ni") )
                .def("get_node_data",[] (AdjList const & ptr, Index ni) -> DataType const & {
                    return ptr.get_node_data(ni); } , py::arg("ni"))
                .def("get_node_data",[] (AdjList  & ptr, Index ni) -> DataType & {
                    return ptr.get_node_data(ni); } , py::arg("ni"))
                .def("get_node",[] (AdjList const & ptr, Index ni) -> Node<DataType> const & {
                    return ptr.get_node(ni); } , py::arg("ni"))
                .def("get_node",[] (AdjList  & ptr, Index ni) -> Node<DataType> & {
                    return ptr.get_node(ni); } , py::arg("ni"))
                .def("get_connected_node_info",[] (AdjList const & ptr, NodeIndexandEdge const & nei) -> NodeIndexandEdge {
                    return ptr.get_connected_node_info(nei); } , py::arg("nei"))
                .def("edge_between_nodes",[] (AdjList const & ptr, Index n1, Index n2) -> bool {
                    return ptr.edge_between_nodes(n1, n2); } , py::arg("n1"), py::arg("n2"))
                .def("edge_index_empty",[] (AdjList const & ptr, Index ni, Index ei) -> bool {
                    return ptr.edge_index_empty(ni, ei); }, py::arg("ni"), py::arg("ei") )
            // operators
//		.def(py::self = py::self)
                ;
        py::class_<UnDirGraph, std::shared_ptr<UnDirGraph>>(m, String(String("UndirectedGraph") + PyType).c_str())
                // ctors
                .def(py::init<>())
                .def(py::init<UnDirGraph const &>(), py::arg("g"))
                ;


        py::class_<Visited, std::shared_ptr<Visited>>(m, String(String("VisitedNode") + PyType).c_str())
                .def(py::init<std::shared_ptr<Visited>, Index>(), py::arg("nparent"), py::arg("nindex") )
                        // methods
                .def("index_in_path",[] (Visited  & ptr, Index i) -> bool {
                    return ptr.index_in_path(i); }, py::arg("i") )
                .def("path_length",[] (Visited  & ptr) -> int {
                    return ptr.path_length(); } )
                        // public attributes
                .def_readwrite("parent", &Visited::parent)
                .def_readwrite("index", &Visited::index)
                ;
        py::class_<Graph, std::shared_ptr<Graph>>(m, String(String("Graph") + PyType).c_str())
                // ctors
                .def(py::init<>())
                        // methods
                .def("begin",[] (Graph  & ptr) { return ptr.begin(); } )
                .def("end",[] (Graph  & ptr) { return ptr.end(); } )
                .def("begin",[] (Graph const & ptr) { return ptr.begin(); } )
                .def("end",[] (Graph const & ptr) { return ptr.end(); } )
                .def("setup_transversal",[] (Graph  & ptr, Index start_n) {
                    ptr.setup_transversal(start_n); }, py::arg("start_n") )
                .def("setup_path_transversal",[] (Graph  & ptr, Index start_n, Index end_n) {
                    ptr.setup_path_transversal(start_n, end_n); }, py::arg("start_n"), py::arg("end_n") )
                .def("add_node",[] (Graph  & ptr, DataType const & d, Size n_edges) ->  Index {
                    return ptr.add_node(d, n_edges); }, py::arg("d"), py::arg("n_edges") )
                .def("add_edge",[] (Graph  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
                    ptr.add_edge(nie1, nie2); }, py::arg("nie1"), py::arg("nie2") )
                .def("remove_node",[] (Graph  & ptr, Index ni) {
                    ptr.remove_node(ni); }, py::arg("ni") )
                .def("remove_edge",[] (Graph  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
                    ptr.remove_edge(nie1, nie2); }, py::arg("nie1"), py::arg("nie2") )
                .def("get_num_nodes",[] (Graph const & ptr) ->  size_t {
                    return ptr.get_num_nodes(); } )
                .def("get_num_edges",[] (Graph const & ptr) ->  size_t {
                    return ptr.get_num_edges(); } )
                .def("get_node_edges",[] (Graph const & ptr, Index ni) ->  std::vector<Edge const *> const & {
                    return ptr.get_node_edges(ni); }, py::arg("ni") )
                .def("get_node",[] (Graph const & ptr, Index ni) ->  Node<DataType> const & {
                    return ptr.get_node(ni); }, py::arg("ni") )
                .def("get_node_data",[] (Graph const & ptr, Index ni) ->  DataType const & {
                    return ptr.get_node_data(ni); }, py::arg("ni") )
                .def("get_node_data",[] (Graph  & ptr, Index ni) ->  DataType & {
                    return ptr.get_node_data(ni); }, py::arg("ni") )
                .def("get_connected_node_info",[] (Graph const & ptr, NodeIndexandEdge const & nei) ->  NodeIndexandEdge {
                    return ptr.get_connected_node_info(nei); }, py::arg("nei") )
                .def("edge_between_nodes",[] (Graph const & ptr, Index n1, Index n2) ->  bool {
                    return ptr.edge_between_nodes(n1, n2); }, py::arg("n1"), py::arg("n2") )
                .def("edge_index_empty",[] (Graph const & ptr, Index ni, Index ei) ->  bool {
                    return ptr.edge_index_empty(ni, ei); }, py::arg("ni"), py::arg("ei") )
                ;
    }
}

namespace data_structure::graph {

    namespace  py = pybind11 ;
    void
    add_bindings(py::module_&);
    template<
            typename DataType,
            typename Graph = Graph<DataType>,
            typename Connection = GraphConnection<DataType>,
            typename Dynamic = GraphDynamic<DataType>,
            typename Iter = GraphIterator<DataType>
            >
    void
    add_graph_specialization(py::module_ & m, String const&  PyType)  {

        // classes
        py::class_<Graph, std::shared_ptr<Graph>>(m, String(String("Graph") + PyType).c_str())
                // ctors
                .def(py::init<>())
                        // methods
                .def("begin",[] (Graph  & ptr)  { return ptr.begin(); } )
                .def("end",[] (Graph  & ptr){ return ptr.end(); } )
                .def("begin",[] (Graph const & ptr) { return ptr.begin(); } )
                .def("end",[] (Graph const & ptr) { return ptr.end(); } )
                .def("transverse",[] (Graph const & ptr, data_structure::graph::GraphNodeOP<DataType> const & graph_node){
                    return ptr.transverse(graph_node); }, py::arg("graph_node") )
                .def("size",[] (Graph const & ptr) ->  size_t { return ptr.size(); } )
                .def("get_node",[] (Graph const & ptr, int index) ->  GraphNodeOP<DataType> const & { return ptr.get_node(index); }, py::arg("index") )
                .def("oldest_node",[] (Graph  & ptr) ->  GraphNodeOP<DataType> { return ptr.oldest_node(); } )
                .def("increase_level",[] (Graph  & ptr) {ptr.increase_level(); } )
                .def("decrease_level",[] (Graph  & ptr) { ptr.decrease_level(); } )
                .def("nodes",[] (Graph const & ptr) ->  GraphNodeOPs<DataType> const & { return ptr.nodes(); } )
                .def("connections",[] (Graph const & ptr) ->  GraphConnectionOPs<DataType> const & { return ptr.connections(); } )
                .def("last_node",[] (Graph  & ptr) ->  GraphNodeOP<DataType> const & { return ptr.last_node(); } )
                .def("level",[] (Graph  & ptr) ->  int { return ptr.level(); } )
                .def("index",[] (Graph  & ptr) ->  int { return ptr.index(); } )
                .def("index",[] (Graph  & ptr, int nindex) { ptr.index(nindex); }, py::arg("nindex") )
                ;

        py::class_<Connection, std::shared_ptr<Connection>>(m, String(String("GraphConnection") + PyType).c_str())
                // ctors
                .def(py::init<GraphNodeOP<DataType>,GraphNodeOP<DataType>,int,int>(),
                     py::arg("node_1"), py::arg("node_2"), py::arg("end_index_1"), py::arg("end_index_2"))
                        // methods
                .def("disconnect",[] (Connection  & ptr) { ptr.disconnect(); } )
                .def("partner",[] (Connection  & ptr, int i) -> GraphNodeOP<DataType> const & {
                    return ptr.partner(i); }, py::arg("i") )
                .def("end_index",[] (Connection  & ptr, int n_index) -> int {
                    return ptr.end_index(n_index); }, py::arg("n_index") )
                .def("node_1",[] (Connection  & ptr) ->  GraphNodeOP<DataType> const & {
                    return ptr.node_1(); } )
                .def("node_2",[] (Connection  & ptr) ->  GraphNodeOP<DataType> const & {
                    return ptr.node_2(); } )
                .def("end_index_1",[] (Connection  & ptr) ->  int {
                    return ptr.end_index_1(); } )
                .def("end_index_2",[] (Connection  & ptr) ->  int {
                    return ptr.end_index_2(); } )
                ;
        py::class_<Dynamic, std::shared_ptr<Dynamic>>(m, String(String("GraphDynamic") + PyType).c_str())
                // ctors
                .def(py::init<>())
                        // methods
                .def("add_data",[] (Dynamic  & ptr, DataType const & data, int parent_index, int orphan) ->  int {
                    return ptr.add_data(data, parent_index, orphan); }, py::arg("data"), py::arg("parent_index"), py::arg("orphan") )
                .def("connect",[] (Dynamic  & ptr, int i, int j) {
                    ptr.connect(i, j); }, py::arg("i"), py::arg("j") )
                ;

        py::class_<Iter, std::shared_ptr<Iter>>(m, String( String("GraphIterator") + PyType ).c_str())
                // ctors
                .def(py::init<>())
                        // operators
//		.def(py::self ++ )
                .def(py::self == py::self)
                .def(py::self != py::self)
                ;


        py::class_<GraphNodeCompare<DataType>, std::shared_ptr<GraphNodeCompare<DataType>>>(m, String(String("GraphNodeCompare") + PyType).c_str())
                .def(py::init())
                        // operators
                .def("compare", [] (GraphNodeCompare<DataType> g, GraphNodeOP<DataType> const& node1,
                                    GraphNodeOP<DataType> const&   node2){
                    return g(node1, node2);
                })
//		.def(py::self () GraphNodeOP<DataType> const &)
                ;

        py::class_<GraphNodeDynamic<DataType>, std::shared_ptr<GraphNodeDynamic<DataType>>>(m, String(String("GraphNodeDynamic") + PyType).c_str())
                // ctors
                .def(py::init<DataType const &,int,int>(),
                     py::arg("data"), py::arg("index"), py::arg("level"))
                        // methods
                .def("add_connection",[] (GraphNodeDynamic<DataType>  & ptr, data_structure::graph::GraphConnectionOP<DataType> const & connection, int pos) {
                    ptr.add_connection(connection, pos); }, py::arg("connection"), py::arg("pos") )
                .def("remove_connection",[] (GraphNodeDynamic<DataType>  & ptr, data_structure::graph::GraphConnectionOP<DataType> const & connection) {
                    ptr.remove_connection(connection); }, py::arg("connection") )
                ;

        py::class_<GraphNodeStatic<DataType>, std::shared_ptr<GraphNodeStatic<DataType>>>(m, String(String("GraphNodeStatic") + PyType).c_str())
                // ctors
                .def(py::init<DataType const &,int,int,int>(),
                     py::arg("data"), py::arg("index"), py::arg("level"), py::arg("n_children"))
                .def(py::init<GraphNode<DataType> const &>(), py::arg("n"))
                        // methods
                .def("add_connection",[] (GraphNodeStatic<DataType>  & ptr, GraphConnectionOP<DataType> const & connection, int pos) {
                    ptr.add_connection(connection, pos); }, py::arg("connection"), py::arg("pos") )
                .def("remove_connection",[] (GraphNodeStatic<DataType>  & ptr, GraphConnectionOP<DataType> const & connection) {
                    ptr.remove_connection(connection); }, py::arg("connection") )
                ;


        py::class_<GraphStatic<DataType>, std::shared_ptr<GraphStatic<DataType>>>(m, String(String("GraphStatic") + PyType).c_str())
                // ctors
                .def(py::init<>())
                .def(py::init<GraphStatic<DataType> const &>(), py::arg("g"))
                        // methods
                .def("add_data",[] (GraphStatic<DataType> & ptr, DataType const & data, int parent_index, int parent_pos, int child_pos, int n_children, int orphan, int index) ->  int {
                         return ptr.add_data(data, parent_index, parent_pos, child_pos, n_children, orphan, index); },
                     py::arg("data") = -1, py::arg("parent_index") = -1, py::arg("parent_pos") = -1, py::arg("child_pos") = -1,
                     py::arg("n_children") = 0, py::arg("orphan") = 0, py::arg("index") = -1
                )
                .def("connect",[] (GraphStatic<DataType> & ptr, int i, int j, int i_pos, int j_pos) {
                    ptr.connect(i, j, i_pos, j_pos); }, py::arg("i"), py::arg("j"), py::arg("i_pos"), py::arg("j_pos") )
                .def("check_pos_is_valid",[] (GraphStatic<DataType>  & ptr, data_structure::graph::GraphNodeOP<DataType> const & n, int & pos) ->  int {
                    return ptr.check_pos_is_valid(n, pos); }, py::arg("n"), py::arg("pos") )
                .def("get_available_pos",[] (GraphStatic<DataType>  & ptr, data_structure::graph::GraphNodeOP<DataType> const & n, int & pos) ->  Ints {
                    return ptr.get_available_pos(n, pos); } , py::arg("n"), py::arg("pos"))
                .def("remove_node",[] (GraphStatic<DataType>  & ptr, int pos) {
                    ptr.remove_node(pos); }, py::arg("pos"))
                .def("remove_level",[] (GraphStatic<DataType>  & ptr, int level) {
                    ptr.remove_level(level); }, py::arg("level") )
                ;

    }
}

namespace data_structure::tree {
    namespace  py = pybind11 ;
    void
    add_bindings(py::module_&);

    template<
            typename DataType,
            typename Tree = Tree<DataType>,
            typename TreeDynamic = TreeDynamic<DataType>,
            typename TreeStatic = TreeStatic<DataType>,
            typename TreeNodeDynamic  = TreeNodeDynamic<DataType>,
            typename TreeNodeStatic = TreeNodeStatic<DataType>
            >
    void
    add_graph_specialization(py::module_ & m, String const& PyType) {

        // classes
        py::class_<Tree, std::shared_ptr<Tree>>(m, String( String("Tree") + PyType).c_str())
                // ctors
                .def(py::init<>())
                        // methods
                .def("begin",[] (Tree  & ptr) { return ptr.begin(); } )
                .def("end",[] (Tree  & ptr) { return ptr.end(); } )
                .def("begin",[] (Tree const & ptr) { return ptr.begin(); } )
                .def("end",[] (Tree const & ptr) { return ptr.end(); } )
                .def("get_node",[] (Tree  & ptr, int index) ->  TreeNodeOP<DataType> const & {
                    return ptr.get_node(index); } , py::arg("index"))
                .def("remove_node",[] (Tree  & ptr, data_structure::tree::TreeNodeOP<DataType> const & n) {
                    ptr.remove_node(n); }, py::arg("n") )
                .def("remove_node",[] (Tree  & ptr, int pos) {
                    ptr.remove_node(pos); }, py::arg("pos") )
                .def("remove_level",[] (Tree  & ptr, int level) {
                    ptr.remove_level(level); }, py::arg("level") )
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

        py::class_<TreeDynamic, std::shared_ptr<TreeDynamic>>(m, String(String("TreeDynamic") + PyType).c_str())
                // ctors
                .def(py::init<>())
                        // methods
                .def("add_data",[] (TreeDynamic  & ptr, DataType const & data, int parent_index) ->  int {
                    return ptr.add_data(data, parent_index); }, py::arg("data"), py::arg("parent_index") = -1 )
                ;

        py::class_<TreeNodeDynamic, std::shared_ptr<TreeNodeDynamic>>(m, String( String("TreeNodeDynamic") + PyType).c_str())
                // ctors
                .def(py::init<DataType const &,int,int>(),
                     py::arg("data"), py::arg("index"), py::arg("level"))
                        // methods
                .def("add_child",[] (TreeNodeDynamic  & ptr, data_structure::tree::TreeNodeOP<DataType> const & c, int pos) {
                    ptr.add_child(c, pos); }, py::arg("c"), py::arg("pos") = -1 )
                .def("remove_child",[] (TreeNodeDynamic  & ptr, data_structure::tree::TreeNodeOP<DataType> const & child) {
                    ptr.remove_child(child); }, py::arg("child") )
                .def("leaf",[] (TreeNodeDynamic  & ptr) ->  bool {
                    return ptr.leaf(); } )
                ;

        py::class_<TreeNodeStatic, std::shared_ptr<TreeNodeStatic>>(m, String( String("TreeNodeStatic") + PyType).c_str() )
                // ctors
                .def(py::init<DataType const &,int,int,int>(),
                     py::arg("data"), py::arg("index"), py::arg("level"), py::arg("n_children"))
                        // methods
                .def("add_child",[] (TreeNodeStatic  & ptr, data_structure::tree::TreeNodeOP<DataType> const & child, int pos) {
                    ptr.add_child(child, pos); }, py::arg("child"), py::arg("pos") = -1 )
                .def("remove_child",[] (TreeNodeStatic  & ptr, data_structure::tree::TreeNodeOP<DataType> const & child) {
                    ptr.remove_child(child); }, py::arg("child") )
                .def("leaf",[] (TreeNodeStatic  & ptr) ->  bool {
                    return ptr.leaf(); } )
                ;

        py::class_<TreeStatic, std::shared_ptr<TreeStatic>>(m, String( String( "TreeStatic") + PyType ).c_str() )
                // ctors
                .def(py::init<>())
                .def(py::init<TreeStatic const &>(), py::arg("t"))
                        // methods
                .def("add_data",[] (TreeStatic  & ptr, DataType const & data, int n_children, int parent_index, int parent_child_index) ->  int {
                         return ptr.add_data(data, n_children, parent_index, parent_child_index); },
                     py::arg("data"), py::arg("n_children") = 1, py::arg("parent_index") = -1, py::arg("parent_child_index") = -1)
                .def("get_available_pos",[] (TreeStatic  & ptr, data_structure::tree::TreeNodeOP<DataType> const & n, int pos) ->  Ints {
                    return ptr.get_available_pos(n, pos); }, py::arg("n"), py::arg("pos") = -1 )
                .def("index",[] (TreeStatic  & ptr, int n_index) {
                    ptr.index(n_index); }, py::arg("n_index") )
                ;

    }
}

namespace data_structure{
    namespace py = pybind11;

    void
    add_bindings(py::module_ & m) {
        auto graph = m.def_submodule("graph");
        graph::add_bindings(graph);
        auto tree = m.def_submodule("tree");
        tree::add_bindings(tree);

        add_graph_specialization<int>(m, "Int");
//        add_graph_specialization<float>(m, "Float");

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// data_structure
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // exceptions
        py::register_exception<GraphException>(m, "GraphException");

        py::class_<DynamicEdges, std::shared_ptr<DynamicEdges>>(m, "DynamicEdges")
                ;

        py::class_<Edge, std::shared_ptr<Edge>>(m, "Edge")
                // ctors
                .def(py::init<Index,Index,Index,Index>(),
                        py::arg("nnode_i"), py::arg("nnode_j"), py::arg("nedge_i"), py::arg("nedge_j"))
                        // methods
                .def("partner",[] (Edge const & ptr, Index index) -> Index {
                    return ptr.partner(index); }, py::arg("index") )
                .def("end_index",[] (Edge const & ptr, Index index) -> Index {
                    return ptr.end_index(index); } , py::arg("index"))
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

        py::class_<NodeIndexandEdgeCompare, std::shared_ptr<NodeIndexandEdgeCompare>>(m, "NodeIndexandEdgeCompare")
	        .def(py::init())
            // operators
	        .def("compare", [](NodeIndexandEdgeCompare& ptr, const NodeIndexandEdge& nie1, const NodeIndexandEdge& nie2) {
	            return ptr(nie1, nie2);
	        }, py::arg("nie1"), py::arg("nie2"))
		;

    }
}

namespace data_structure::graph {
    namespace py = pybind11;
    void
    add_bindings(py::module_& m) {
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// data_structure::graph
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        py::register_exception<GraphException>(m, "GraphException");
        py::enum_<GraphNodeType>(m, "GraphNodeType")
            .value("GraphNodeTypeStatic", GraphNodeType::GraphNodeTypeStatic)
            .value("GraphNodeTypeDynamic", GraphNodeType::GraphNodeTypeDynamic) ;

        add_graph_specialization<int>(m, "Int");
//        add_graph_specialization<float>(m, "Float");
    }
}

namespace data_structure::tree {
    namespace py = pybind11;
    void
    add_bindings(py::module_ & m) {
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// data_structure::tree
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        // Errors
        py::register_exception<TreeException>(m, "TreeException");

        add_graph_specialization<int>(m, "Int");
        add_graph_specialization<float>(m, "Float");

    }
}

#endif // PYBIND11_DATA_STRUCTURE_HPP