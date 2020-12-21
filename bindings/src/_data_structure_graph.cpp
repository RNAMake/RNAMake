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

// motif includes
#include <motif/motif.h>
#include <motif/motif_ensemble.h>
#include <motif/motif_factory.h>
#include <motif/motif_scorer.h>
#include <motif/motif_state.h>
#include <motif/motif_state_aligner.h>
#include <motif/motif_state_ensemble.h>
#include <motif/motif_to_secondary_structure.h>
#include <motif/pose.h>
#include <motif/pose_factory.h>

// motif_data_structure includes
#include <motif_data_structure/motif_connection.h>
#include <motif_data_structure/motif_graph.h>
#include <motif_data_structure/motif_merger.h>
#include <motif_data_structure/motif_state_ensemble_graph.h>
#include <motif_data_structure/motif_state_ensemble_tree.h>
#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_data_structure/motif_state_node.hpp>
#include <motif_data_structure/motif_state_tree.h>
#include <motif_data_structure/motif_topology.h>
#include <motif_data_structure/motif_tree.h>

// motif_search includes
#include <motif_search/motif_state_monte_carlo.h>
#include <motif_search/problem.h>
#include <motif_search/search.h>
#include <motif_search/solution_filter.h>
#include <motif_search/solution_topology.h>

// motif_search::exhaustive includes
#include <motif_search/exhaustive/motif_state_enumerator.h>
#include <motif_search/exhaustive/scorer.h>
#include <motif_search/exhaustive/search.h>

// motif_search::monte_carlo includes
#include <motif_search/monte_carlo/scorer.h>
#include <motif_search/monte_carlo/search.h>

// motif_search::path_finding includes
#include <motif_search/path_finding/node.h>
#include <motif_search/path_finding/scorer.h>
#include <motif_search/path_finding/search.h>
#include <motif_search/path_finding/selector.h>

// motif_tools includes
#include <motif_tools/segmenter.h>

// resources includes
#include <resources/added_motif_library.h>
#include <resources/motif_ensemble_sqlite_connection.h>
#include <resources/motif_sqlite_connection.h>
#include <resources/motif_sqlite_library.h>
#include <resources/motif_state_ensemble_sqlite_library.h>
#include <resources/motif_state_sqlite_library.h>
#include <resources/resource_manager.h>
#include <resources/sqlite_library.h>

// secondary_structure includes
#include <secondary_structure/basepair.h>
#include <secondary_structure/chain.h>
#include <secondary_structure/motif.h>
#include <secondary_structure/pose.h>
#include <secondary_structure/residue.h>
#include <secondary_structure/rna_structure.h>
#include <secondary_structure/secondary_structure_parser.h>
#include <secondary_structure/secondary_structure_tree.h>
#include <secondary_structure/sequence_constraint.h>
#include <secondary_structure/sequence_tools.h>
#include <secondary_structure/structure.h>
#include <secondary_structure/util.h>

// sequence_optimization includes
#include <sequence_optimization/sequence_optimizer_3d.hpp>

// structure includes
#include <structure/atom.h>
#include <structure/basepair.h>
#include <structure/basepair_state.h>
#include <structure/beads.h>
#include <structure/chain.h>
#include <structure/cif_parser.h>
#include <structure/close_chain.h>
#include <structure/is_equal.h>
#include <structure/pdb_parser.h>
#include <structure/residue.h>
#include <structure/residue_type.h>
#include <structure/residue_type_set.h>
#include <structure/residue_type_set_manager.h>
#include <structure/rna_structure.h>
#include <structure/structure.h>

// thermo_fluctuation includes
#include <thermo_fluctuation/thermo_fluc_relax.h>
#include <thermo_fluctuation/thermo_fluc_sampler.h>
#include <thermo_fluctuation/thermo_fluc_scorer.h>
#include <thermo_fluctuation/thermo_fluc_simulation.h>

// thermo_fluctuation::graph includes
#include <thermo_fluctuation/graph/sampler.h>
#include <thermo_fluctuation/graph/scorer.h>
#include <thermo_fluctuation/graph/simulation.h>

// thermo_fluctuation::graph::logging includes
#include <thermo_fluctuation/graph/logging.h>

// thermo_fluctuation::graph::sterics includes
#include <thermo_fluctuation/graph/sterics.h>

// util includes
#include <util/basic_io.hpp>
#include <util/cartesian_product.h>
#include <util/monte_carlo.h>
#include <util/motif_type.h>
#include <util/random_number_generator.h>
#include <util/sqlite3_connection.h>
#include <util/steric_lookup.hpp>
#include <util/uuid.h>
#include <util/x3dna.h>
using namespace data_structure::graph;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_data_structure_graph,m) {

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
}