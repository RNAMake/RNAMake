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
using namespace data_structure::tree;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_data_structure_tree,m) {

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