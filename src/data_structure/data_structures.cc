//
// Created by Joseph Yesselman on 10/22/17.
//

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

//rnamake headers
#include <base/types.h>
#include <data_structure/graph.h>

#include <primitives/residue.h>
#include <primitives/chain.h>

namespace py = pybind11;
namespace data_structure {

PYBIND11_PLUGIN(data_structure) {
    py::module m("data_structure", "basic organization classes");
    // NodeIndexandEdge
    py::class_<NodeIndexandEdge, std::shared_ptr<NodeIndexandEdge>>(m, "NodeIndexandEdge")
            .def(py::init<Index, Index>())
            .def_readonly("node_index", &NodeIndexandEdge::node_index)
            .def_readonly("edge_index", &NodeIndexandEdge::edge_index);

    // Edge Class
    py::class_<Edge, std::shared_ptr<Edge>>(m, "Edge")
            .def(py::init<Index, Index, Index, Index>())
            .def_readonly("node_i", &Edge::node_i)
            .def_readonly("node_j", &Edge::node_j)
            .def_readonly("edge_i", &Edge::edge_i)
            .def_readonly("edge_j", &Edge::edge_j);

    // Node Class
    typedef Node<int> NodeInt;
    py::class_<NodeInt, std::shared_ptr<NodeInt>>(m, "NodeInt")
            .def(py::init<int, Index>())
            .def("data", (int & (NodeInt::*)(void)) &NodeInt::data)
            .def("index", &NodeInt::index);

    // Int Graphs
    typedef FixedEdgeDirectedGraph<int> FixedEdgeDirectedGraphInt;
    py::class_<FixedEdgeDirectedGraphInt, std::shared_ptr<FixedEdgeDirectedGraphInt> >(m, "FixedEdgeDirectedGraphInt")
            .def(py::init<>())
            .def("__len__", &FixedEdgeDirectedGraphInt::get_num_nodes)
            .def("__iter__", [](FixedEdgeDirectedGraphInt const & g) {
                    return py::make_iterator(g.begin(), g.end()); }, py::keep_alive<0, 1>())
            .def("setup_transversal", &FixedEdgeDirectedGraphInt::setup_transversal)
            .def("setup_path_transversal", &FixedEdgeDirectedGraphInt::setup_path_transversal)
            .def("setup_sub_graph_transversal", &FixedEdgeDirectedGraphInt::setup_sub_graph_transversal)
            .def("add_node", (int (FixedEdgeDirectedGraphInt::*)(int const &, Size))
                    &FixedEdgeDirectedGraphInt::add_node)
            .def("add_node", (int (FixedEdgeDirectedGraphInt::*)(int const &, Size, Index, NodeIndexandEdge const &)) &FixedEdgeDirectedGraphInt::add_node)
            .def("add_edge", &FixedEdgeDirectedGraphInt::add_edge)
            .def("remove_node", &FixedEdgeDirectedGraphInt::remove_node)
            .def("remove_edge", &FixedEdgeDirectedGraphInt::remove_edge)
            .def("get_num_nodes", &FixedEdgeDirectedGraphInt::get_num_nodes)
            .def("get_num_edges", &FixedEdgeDirectedGraphInt::get_num_edges)
            .def("get_node_edges", &FixedEdgeDirectedGraphInt::get_num_edges)
            .def("get_node", &FixedEdgeDirectedGraphInt::get_node)
            .def("get_node_data", (int & (FixedEdgeDirectedGraphInt::*)(int)) &FixedEdgeDirectedGraphInt::get_node_data)
            .def("get_connected_node_info", &FixedEdgeDirectedGraphInt::get_connected_node_info)
            .def("edge_between_nodes", &FixedEdgeDirectedGraphInt::edge_between_nodes)
            .def("edge_index_empty", &FixedEdgeDirectedGraphInt::edge_index_empty);


    // Graph Class
    //typedef primitives::Residue DataType;
    //typedef std::shared_ptr<DataType> DataTypeOP;
    //typedef Graph<DataType> ResidueGraph;
    /*py::class_<ResidueGraph, std::shared_ptr<ResidueGraph> >(m, "ResidueGraph")
            .def(py::init<>())
            .def("__len__", &ResidueGraph::get_num_nodes)
            .def("__iter__", [](ResidueGraph const & g) { return py::make_iterator(g.begin(), g.end()); },
                 py::keep_alive<0, 1>())
            .def("change_transverse_start", &ResidueGraph::change_transverse_start)
            .def("get_edges", &ResidueGraph::get_edges)
            .def("get_num_edges", &ResidueGraph::get_num_edges)
            .def("get_num_nodes", &ResidueGraph::get_num_nodes)
            .def("are_nodes_connected", &ResidueGraph::are_nodes_connected)
            .def("get_node", &ResidueGraph::get_node)
            .def("add_node", (int (ResidueGraph::*)(DataTypeOP)) &ResidueGraph::add_node)
            .def("add_node", (int (ResidueGraph::*)(DataTypeOP, Index)) &ResidueGraph::add_node);

    py::class_<ResidueGraph::Node, std::shared_ptr<ResidueGraph::Node> >(m, "ResidueGraphNode")
            .def(py::init<Index, DataType const *>())
            .def_readonly("index", &ResidueGraph::Node::index)
            .def_readonly("data", &ResidueGraph::Node::data);

    typedef base::VectorContainer<Edge const *> Edges;
    typedef std::shared_ptr<Edges> EdgesOP;
    py::class_<Edges, EdgesOP>(m, "Basepairs")
            .def(py::init<std::vector<Edge const *> const &>())
    .def("__iter__", []( Edges const & c) {
        return py::make_iterator(c.begin(), c.end()); }, py::keep_alive<0, 1>())
            .def("__len__", &Edges::size)
            .def("size", &Edges::size)
            .def("get_data", &Edges::get_data)
            .def("__getitem__", &Edges::operator[]);

    typedef DirectedGraph<DataType> ResidueDirectedGraph;
    py::class_<ResidueDirectedGraph, std::shared_ptr<ResidueDirectedGraph> >(m, "ResidueDirectedGraph")
            .def(py::init<>())
            .def("__len__", &ResidueDirectedGraph::get_num_nodes)
            .def("__iter__", [](ResidueDirectedGraph const & g) { return py::make_iterator(g.begin(), g.end()); },
                 py::keep_alive<0, 1>())
            .def("change_transverse_start", &ResidueDirectedGraph::change_transverse_start)
            .def("get_edges", &ResidueDirectedGraph::get_edges)
            .def("get_num_edges", &ResidueDirectedGraph::get_num_edges)
            .def("get_num_nodes", &ResidueDirectedGraph::get_num_nodes)
            .def("are_nodes_connected", &ResidueDirectedGraph::are_nodes_connected)
            .def("get_node", &ResidueDirectedGraph::get_node)
            .def("add_node", (int (ResidueDirectedGraph::*)(DataTypeOP, Size)) &ResidueDirectedGraph::add_node)
            .def("add_node", (int (ResidueDirectedGraph::*)(DataTypeOP, Size, Index, Index, Index)) &ResidueDirectedGraph::add_node)
            .def("edge_index_empty", &ResidueDirectedGraph::edge_index_empty)
            .def("add_edge", &ResidueDirectedGraph::add_edge)
            .def("get_roots", &ResidueDirectedGraph::get_roots)
            .def("get_parent_index", &ResidueDirectedGraph::get_parent_index)
            .def("has_parent", &ResidueDirectedGraph::has_parent)
            .def("get_parent_edge_index", &ResidueDirectedGraph::get_parent_edge_index)
            .def("get_all_edges", &ResidueDirectedGraph::get_all_edges)
            .def("edge_exists", &ResidueDirectedGraph::edge_exist);*/


    /*typedef primitives::PrimitiveChain DataType1;
    typedef std::shared_ptr<DataType1> DataTypeOP1;
    typedef Graph<DataType1> GraphType1;
    py::class_<GraphType1, std::shared_ptr<GraphType1> >(m, "Graph")
            .def(py::init<>())
            .def("__len__", &GraphType1::get_num_nodes);*/

    // StaticEdgedGraph Class
    /*typedef StaticEdgedGraph<DataType> StaticEdgedGraphType;
    py::class_<StaticEdgedGraphType, std::shared_ptr<StaticEdgedGraphType> >(m, "StaticEdgedGraph")
            .def(py::init<>())
            .def("__len__", &StaticEdgedGraphType::get_num_nodes)
            .def("__iter__", [](StaticEdgedGraphType const & g) { return py::make_iterator(g.begin(), g.end()); },
                 py::keep_alive<0, 1>())
            //.def("change_transverse_start", &StaticEdgedGraphType::change_transverse_start)
            .def("get_edges", &StaticEdgedGraphType::get_edges)
            .def("get_num_edges", &StaticEdgedGraphType::get_num_edges)
            .def("get_num_nodes", &StaticEdgedGraphType::get_num_nodes)
            .def("are_nodes_connected", &StaticEdgedGraphType::are_nodes_connected)
            .def("get_node", &StaticEdgedGraphType::get_node)
            .def("add_node", (int (StaticEdgedGraphType::*)(DataTypeOP)) &StaticEdgedGraphType::add_node)
            .def("add_node", (int (StaticEdgedGraphType::*)(DataTypeOP, Index)) &StaticEdgedGraphType::add_node);
    */

    return m.ptr();

}
}

