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
using namespace motif_data_structure;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_motif_data_structure,m) {

	// free functions

    m.def("graph_to_tree", [] (MotifGraphOP const & mg, data_structure::graph::GraphNodeOP<motif::MotifOP> start, structure::BasepairOP last_end) -> MotifTreeOP {
		 return graph_to_tree(mg, start, last_end); },
	py::arg("mg"),
	py::arg("start") = nullptr,
	py::arg("last_end") = nullptr 
    );

    m.def("motif_state_ensemble_graph_from_motif_graph", [] (MotifGraph & mg, resources::Manager & man, MotifStateEnsembleGraph & meg, std::map<int, int> & map) {
		motif_state_ensemble_graph_from_motif_graph(mg, man, meg, map); },
	py::arg("graph"),
	py::arg("manager"),
	py::arg("ensemble_graph"),
	py::arg("map")
    );
	// classes

        py::class_<ChainNodeData, std::shared_ptr<ChainNodeData>>(m, "ChainNodeData")
		// ctors
		.def(py::init<>())
		.def(py::init<structure::ChainOP const &,util::Uuid const &,int,int>())
		.def(py::init<ChainNodeData const &>())
		// methods
		.def("included_res",[] (ChainNodeData  & ptr) -> structure::ResidueOPs {
		 return ptr.included_res(); } )
		// public attributes
		.def_readwrite("c", &ChainNodeData::c)
		.def_readwrite("m_id", &ChainNodeData::m_id)
		.def_readwrite("prime5_override", &ChainNodeData::prime5_override)
		.def_readwrite("prime3_override", &ChainNodeData::prime3_override)
		;

        py::class_<GraphtoTree, std::shared_ptr<GraphtoTree>>(m, "GraphtoTree")
		// ctors
		.def(py::init<>())
		// methods
		.def("convert",[] (GraphtoTree  & ptr, MotifGraphOP const & mg, data_structure::graph::GraphNodeOP<motif::MotifOP> start, int start_end_index, data_structure::graph::GraphNodeOP<motif::MotifOP> last_node) -> MotifTreeOP {
		 return ptr.convert(mg, start, start_end_index, last_node); } )
		;

        py::class_<MSNodeData, std::shared_ptr<MSNodeData>>(m, "MSNodeData")
		// ctors
		.def(py::init<motif::MotifStateOP const &>())
		.def(py::init<MSNodeData const &>())
		// methods
		.def("get_end_state",[] (MSNodeData  & ptr, String const & name) ->  structure::BasepairStateOP {
		 return ptr.get_end_state(name); } )
		.def("get_end_state",[] (MSNodeData  & ptr, int i) ->  structure::BasepairStateOP {
		 return ptr.get_end_state(i); } )
		.def("get_end_index",[] (MSNodeData  & ptr, String const & name) ->  int {
		 return ptr.get_end_index(name); } )
		.def("name",[] (MSNodeData  & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("block_end_add",[] (MSNodeData  & ptr) ->  int {
		 return ptr.block_end_add(); } )
		.def("end_name",[] (MSNodeData  & ptr, int i) ->  String const & {
		 return ptr.end_name(i); } )
		.def("uuid",[] (MSNodeData  & ptr) ->  util::Uuid const & {
		 return ptr.uuid(); } )
		.def("uuid",[] (MSNodeData  & ptr, util::Uuid const & uuid) {
		ptr.uuid(uuid); } )
		// public attributes
		.def_readwrite("ref_state", &MSNodeData::ref_state)
		.def_readwrite("cur_state", &MSNodeData::cur_state)
		;

        py::class_<MotifConnection, std::shared_ptr<MotifConnection>>(m, "MotifConnection")
		// ctors
		.def(py::init<>())
		.def(py::init<int,int,String const &,String const &>())
		.def(py::init<MotifConnection const &>())
		// methods
		.def("i",[] (MotifConnection  & ptr) ->  int {
		 return ptr.i(); } )
		.def("j",[] (MotifConnection  & ptr) ->  int {
		 return ptr.j(); } )
		.def("name_i",[] (MotifConnection  & ptr) ->  String const & {
		 return ptr.name_i(); } )
		.def("name_j",[] (MotifConnection  & ptr) ->  String const & {
		 return ptr.name_j(); } )
		.def("name_i",[] (MotifConnection  & ptr, String const & name) {
		ptr.name_i(name); } )
		.def("name_j",[] (MotifConnection  & ptr, String const & name) {
		ptr.name_j(name); } )
		.def("to_str",[] (MotifConnection  & ptr) -> String {
		 return ptr.to_str(); } )
		;

        py::class_<MotifConnections, std::shared_ptr<MotifConnections>>(m, "MotifConnections")
		// ctors
		.def(py::init<>())
		.def(py::init<MotifConnections const>())
		// methods
		.def("begin",[] (MotifConnections  & ptr) -> MotifConnectionOPs::iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifConnections  & ptr) -> MotifConnectionOPs::iterator {
		 return ptr.end(); } )
		.def("begin",[] (MotifConnections const & ptr) -> MotifConnectionOPs::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifConnections const & ptr) -> MotifConnectionOPs::const_iterator {
		 return ptr.end(); } )
		.def("size",[] (MotifConnections const & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("add_connection",[] (MotifConnections  & ptr, int i, int j, String const & name_i, String const & name_j) {
		ptr.add_connection(i, j, name_i, name_j); } )
		.def("remove_connections_to",[] (MotifConnections  & ptr, int index) {
		ptr.remove_connections_to(index); } )
		.def("in_connection",[] (MotifConnections const & ptr, int index, String const & name) -> bool {
		 return ptr.in_connection(index, name); } )
		.def("update_connection_name",[] (MotifConnections  & ptr, int index, String const & name, String const & new_name) {
		ptr.update_connection_name(index, name, new_name); } )
		;

        py::class_<MotifGraph, std::shared_ptr<MotifGraph>>(m, "MotifGraph")
		// ctors
		.def(py::init<>())
		//.def(py::init<String const>())
		.def(py::init<String const,MotifGraphStringType const>())
		.def(py::init<MotifGraph const>())
		// methods
		.def("update_indexes",[] (MotifGraph  & ptr, std::map<int, int> const & index_hash) {
		ptr.update_indexes(index_hash); } )
		.def("begin",[] (MotifGraph  & ptr) -> MotifGraph::iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifGraph  & ptr) -> MotifGraph::iterator {
		 return ptr.end(); } )
		.def("add_motif",[] (MotifGraph  & ptr, motif::MotifOP const & m, int parent_index, int parent_end_index, int orphan) -> int {
		 return ptr.add_motif(m, parent_index, parent_end_index, orphan); } )
		.def("add_motif",[] (MotifGraph  & ptr, motif::MotifOP const & m, int parent_index, String const & p_end_name) -> int {
		 return ptr.add_motif(m, parent_index, p_end_name); } )
		.def("add_motif_tree",[] (MotifGraph  & ptr, MotifTreeOP const & mt, int parent_index, int parent_end_index) {
		ptr.add_motif_tree(mt, parent_index, parent_end_index); } )
		.def("add_motif_tree",[] (MotifGraph  & ptr, MotifTreeOP const & mt, int parent_index, String const & parent_end_name) {
		ptr.add_motif_tree(mt, parent_index, parent_end_name); } )
		.def("add_motif_graph",[] (MotifGraph  & ptr, MotifGraph & mg, int parent_index, int parent_end_index) {
		ptr.add_motif_graph(mg, parent_index, parent_end_index); } )
		.def("add_connection",[] (MotifGraph  & ptr, int i, int j, String const & i_bp_name, String const & j_bp_name) {
		ptr.add_connection(i, j, i_bp_name, j_bp_name); } )
		.def("parent_index",[] (MotifGraph  & ptr, int node_index) ->  int {
		 return ptr.parent_index(node_index); } )
		.def("parent_end_index",[] (MotifGraph  & ptr, int node_index) ->  int {
		 return ptr.parent_end_index(node_index); } )
		.def("remove_motif",[] (MotifGraph  & ptr, int index) {
		ptr.remove_motif(index); } )
		.def("remove_level",[] (MotifGraph  & ptr, int level) {
		ptr.remove_level(level); } )
		.def("size",[] (MotifGraph  & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("increase_level",[] (MotifGraph  & ptr) {
		ptr.increase_level(); } )
		.def("last_node",[] (MotifGraph  & ptr) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> const & {
		 return ptr.last_node(); } )
		.def("oldest_node",[] (MotifGraph  & ptr) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> {
		 return ptr.oldest_node(); } )
		.def("get_node",[] (MotifGraph const & ptr, int i) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> const & {
		 return ptr.get_node(i); } )
		.def("get_node",[] (MotifGraph const & ptr, util::Uuid const & uuid) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> const {
		 return ptr.get_node(uuid); } )
		.def("get_node",[] (MotifGraph const & ptr, String const & m_name) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> const {
		 return ptr.get_node(m_name); } )
		.def("set_index",[] (MotifGraph  & ptr, int index) {
		ptr.set_index(index); } )
		.def("replace_ideal_helices",[] (MotifGraph  & ptr) {
		ptr.replace_ideal_helices(); } )
		.def("replace_ideal_helices",[] (MotifGraph  & ptr) {
		ptr.replace_ideal_helices(); } )
		.def("replace_helical_sequence",[] (MotifGraph  & ptr, secondary_structure::PoseOP const & ss) {
		ptr.replace_helical_sequence(ss); } )
		.def("replace_helical_sequence",[] (MotifGraph  & ptr, String const & seq) {
		ptr.replace_helical_sequence(seq); } )
		.def("designable_secondary_structure",[] (MotifGraph  & ptr) -> secondary_structure::PoseOP {
		 return ptr.designable_secondary_structure(); } )
		.def("designable_sequence",[] (MotifGraph  & ptr) ->  String {
		 return ptr.designable_sequence(); } )
		.def("write_pdbs",[] (MotifGraph  & ptr, String const & fname) {
		ptr.write_pdbs(fname); } )
		.def("topology_to_str",[] (MotifGraph  & ptr) -> String {
		 return ptr.topology_to_str(); } )
		.def("to_str",[] (MotifGraph  & ptr) -> String {
		 return ptr.to_str(); } )
		.def("_update_align_list",[] (MotifGraph  & ptr) {
		ptr._update_align_list(); } )
		.def("_update_merger",[] (MotifGraph  & ptr) {
		ptr._update_merger(); } )
		.def("beads",[] (MotifGraph  & ptr) ->  structure::Beads {
		 return ptr.beads(); } )
		.def("get_available_end",[] (MotifGraph  & ptr, int end) -> structure::BasepairOP const & {
		 return ptr.get_available_end(end); } )
		.def("get_available_end",[] (MotifGraph  & ptr, int pos) -> structure::BasepairOP const & {
		 return ptr.get_available_end(pos); } )
		.def("get_available_end",[] (MotifGraph  & ptr, int pos) -> structure::BasepairOP const &  {
		 return ptr.get_available_end(pos); } )
		.def("unaligned_nodes",[] (MotifGraph const & ptr) -> data_structure::graph::GraphNodeOPs<motif::MotifOP> const {
		 return ptr.unaligned_nodes(); } )
		.def("aligned",[] (MotifGraph  & ptr) -> std::map<int, int> {
		 return ptr.aligned(); } )
		.def("connections",[] (MotifGraph  & ptr) -> data_structure::graph::GraphConnectionOPs<motif::MotifOP> const & {
		 return ptr.connections(); } )
		.def("get_structure",[] (MotifGraph  & ptr) ->  structure::RNAStructureOP const & {
		 return ptr.get_structure(); } )
		.def("secondary_structure",[] (MotifGraph  & ptr) -> secondary_structure::PoseOP {
		 return ptr.secondary_structure(); } )
		.def("to_pdb",[] (MotifGraph  & ptr, String const fname, int renumber, int close_chains, int conect_statement) {
		ptr.to_pdb(fname, renumber, close_chains, conect_statement); } )
		.def("pdb_str",[] (MotifGraph  & ptr, int renumber, int close_chains, int conect_statement) -> String {
		 return ptr.pdb_str(renumber, close_chains, conect_statement); } )
		.def("sequence",[] (MotifGraph  & ptr) ->  String {
		 return ptr.sequence(); } )
		.def("dot_bracket",[] (MotifGraph  & ptr) ->  String {
		 return ptr.dot_bracket(); } )
		.def("get_int_option",[] (MotifGraph  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (MotifGraph  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (MotifGraph  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (MotifGraph  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		.def("set_option_value",[] (MotifGraph  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (MotifGraph  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (MotifGraph  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (MotifGraph  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
		;

        py::class_<MotifMerger, std::shared_ptr<MotifMerger>>(m, "MotifMerger")
		// ctors
		.def(py::init<>())
		.def(py::init<MotifMerger const &,motif::MotifOPs const &>())
		// methods
		.def("add_motif",[] (MotifMerger  & ptr, motif::MotifOP const & m) {
		ptr.add_motif(m); } )
		.def("add_motif",[] (MotifMerger  & ptr, motif::MotifOP const & m) {
		ptr.add_motif(m); } )
		.def("get_structure",[] (MotifMerger  & ptr) -> structure::RNAStructureOP const & {
		 return ptr.get_structure(); } )
		.def("remove_motif",[] (MotifMerger  & ptr, motif::MotifOP const & m) {
		ptr.remove_motif(m); } )
		.def("update_motif",[] (MotifMerger  & ptr, motif::MotifOP const & m) {
		ptr.update_motif(m); } )
		.def("connect_motifs",[] (MotifMerger  & ptr, motif::MotifOP const & m1, motif::MotifOP const & m2, structure::BasepairOP const & m1_end, structure::BasepairOP const & m2_end) {
		ptr.connect_motifs(m1, m2, m1_end, m2_end); } )
		.def("secondary_structure",[] (MotifMerger  & ptr) -> secondary_structure::PoseOP {
		 return ptr.secondary_structure(); } )
		.def("get_residue",[] (MotifMerger  & ptr, util::Uuid const & uuid) ->  structure::ResidueOP {
		 return ptr.get_residue(uuid); } )
		.def("get_basepair",[] (MotifMerger  & ptr, util::Uuid const & uuid) ->  structure::BasepairOP {
		 return ptr.get_basepair(uuid); } )
		.def("to_pdb",[] (MotifMerger  & ptr, String const & fname, int renumber, int close_chains, int conect_statements) {
		ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
		;

        py::class_<MotifStateEnsembleGraph, std::shared_ptr<MotifStateEnsembleGraph>>(m, "MotifStateEnsembleGraph")
		// ctors
		.def(py::init<>())
		.def(py::init<MotifStateEnsembleGraph const &>())
		// methods
		.def("begin",[] (MotifStateEnsembleGraph  & ptr) -> MotifStateEnsembleGraph::iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateEnsembleGraph  & ptr) -> MotifStateEnsembleGraph::iterator {
		 return ptr.end(); } )
		.def("begin",[] (MotifStateEnsembleGraph const & ptr) -> MotifStateEnsembleGraph::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateEnsembleGraph const & ptr) -> MotifStateEnsembleGraph::const_iterator {
		 return ptr.end(); } )
		.def("size",[] (MotifStateEnsembleGraph  & ptr) -> size_t {
		 return ptr.size(); } )
		.def("add_ensemble",[] (MotifStateEnsembleGraph  & ptr, motif::MotifStateEnsemble const & ensemble) -> int {
		 return ptr.add_ensemble(ensemble); } )
		.def("add_ensemble",[] (MotifStateEnsembleGraph  & ptr, motif::MotifStateEnsemble const & ensemble, data_structure::NodeIndexandEdge const & parent_nie) -> int {
		 return ptr.add_ensemble(ensemble, parent_nie); } )
//		.def("add_connection",[] (MotifStateEnsembleGraph  & ptr) -> int {
//		 return ptr.add_connection(); } )
		.def("get_ensemble",[] (MotifStateEnsembleGraph  & ptr, Index ni) ->  motif::MotifStateEnsemble const & {
		 return ptr.get_ensemble(ni); } )
		.def("has_parent",[] (MotifStateEnsembleGraph const & ptr, Index ni) ->  bool {
		 return ptr.has_parent(ni); } )
		.def("get_parent_index",[] (MotifStateEnsembleGraph const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_index(ni); } )
		.def("get_parent_end_index",[] (MotifStateEnsembleGraph const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_end_index(ni); } )
		.def("get_ensemble_connections",[] (MotifStateEnsembleGraph const & ptr, Index ni) ->  std::vector<data_structure::Edge const *> const & {
		 return ptr.get_ensemble_connections(ni); } )
		.def("are_ensembles_connected",[] (MotifStateEnsembleGraph  & ptr, Index n1, Index n2) ->  bool {
		 return ptr.are_ensembles_connected(n1, n2); } )
		;

        py::class_<MotifStateEnsembleOPGraph, std::shared_ptr<MotifStateEnsembleOPGraph>>(m, "MotifStateEnsembleOPGraph")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (MotifStateEnsembleOPGraph  & ptr) -> MotifStateEnsembleOPGraph::iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateEnsembleOPGraph  & ptr) -> MotifStateEnsembleOPGraph::iterator {
		 return ptr.end(); } )
		.def("begin",[] (MotifStateEnsembleOPGraph const & ptr) -> MotifStateEnsembleOPGraph::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateEnsembleOPGraph const & ptr) -> MotifStateEnsembleOPGraph::const_iterator {
		 return ptr.end(); } )
		.def("size",[] (MotifStateEnsembleOPGraph  & ptr) -> size_t {
		 return ptr.size(); } )
		.def("add_ensemble",[] (MotifStateEnsembleOPGraph  & ptr, motif::MotifStateEnsembleOP ensemble) -> int {
		 return ptr.add_ensemble(ensemble); } )
		.def("add_ensemble",[] (MotifStateEnsembleOPGraph  & ptr, motif::MotifStateEnsembleOP ensemble, data_structure::NodeIndexandEdge const & parent_nie) -> int {
		 return ptr.add_ensemble(ensemble, parent_nie); } )
		.def("get_ensemble",[] (MotifStateEnsembleOPGraph  & ptr, Index ni) ->  motif::MotifStateEnsembleOP {
		 return ptr.get_ensemble(ni); } )
		.def("has_parent",[] (MotifStateEnsembleOPGraph const & ptr, Index ni) ->  bool {
		 return ptr.has_parent(ni); } )
		.def("get_parent_index",[] (MotifStateEnsembleOPGraph const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_index(ni); } )
		.def("get_parent_end_index",[] (MotifStateEnsembleOPGraph const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_end_index(ni); } )
		.def("get_leafs",[] (MotifStateEnsembleOPGraph  & ptr) ->  std::vector<data_structure::NodeIndexandEdge> {
		 return ptr.get_leafs(); } )
		;

        py::class_<MotifStateEnsembleTree, std::shared_ptr<MotifStateEnsembleTree>>(m, "MotifStateEnsembleTree")
		// ctors
		.def(py::init<>())
		.def(py::init<MotifTreeOP const>())
		.def(py::init<MotifStateTreeOP const>())
		// methods
		.def("begin",[] (MotifStateEnsembleTree  & ptr) -> MotifStateEnsembleTree::iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateEnsembleTree  & ptr) -> MotifStateEnsembleTree::iterator {
		 return ptr.end(); } )
		.def("begin",[] (MotifStateEnsembleTree const & ptr) -> MotifStateEnsembleTree::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateEnsembleTree const & ptr) -> MotifStateEnsembleTree::const_iterator {
		 return ptr.end(); } )
		.def("add_ensemble",[] (MotifStateEnsembleTree  & ptr, motif::MotifStateEnsembleOP const & ensemble, int parent_index, int parent_end_index) -> int {
		 return ptr.add_ensemble(ensemble, parent_index, parent_end_index); } )
		.def("to_mst",[] (MotifStateEnsembleTree  & ptr) -> MotifStateTreeOP {
		 return ptr.to_mst(); } )
		.def("size",[] (MotifStateEnsembleTree  & ptr) -> size_t {
		 return ptr.size(); } )
		.def("get_node",[] (MotifStateEnsembleTree  & ptr, int i) -> MotifStateEnsembleTreeNodeOP const & {
		 return ptr.get_node(i); } )
		.def("last_node",[] (MotifStateEnsembleTree  & ptr) ->  MotifStateEnsembleTreeNodeOP const & {
		 return ptr.last_node(); } )
		;

        py::class_<MotifStateEnsembleTreeEnumerator, std::shared_ptr<MotifStateEnsembleTreeEnumerator>>(m, "MotifStateEnsembleTreeEnumerator")
		// ctors
		.def(py::init<motif_data_structure::MotifStateEnsembleTreeOP>())
		// methods
		.def("record",[] (MotifStateEnsembleTreeEnumerator  & ptr, String fname) {
		ptr.record(fname); } )
		// public attributes
		.def_readwrite("mtst_", &MotifStateEnsembleTreeEnumerator::mtst_)
		;

        py::class_<MotifStateGraph, std::shared_ptr<MotifStateGraph>>(m, "MotifStateGraph")
		// ctors
		.def(py::init<>())
		.def(py::init<MotifGraphOP const>())
		.def(py::init<MotifStateGraph const>())
		// methods
		.def("begin",[] (MotifStateGraph  & ptr) -> MotifStateGraph::iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateGraph  & ptr) -> MotifStateGraph::iterator {
		 return ptr.end(); } )
		.def("size",[] (MotifStateGraph  & ptr) -> size_t {
		 return ptr.size(); } )
		.def("add_state",[] (MotifStateGraph  & ptr, motif::MotifStateOP const & state, int parent_index, int parent_end_index, int orphan) -> int {
		 return ptr.add_state(state, parent_index, parent_end_index, orphan); } )
		.def("add_state",[] (MotifStateGraph  & ptr, motif::MotifStateOP const & state, int parent_index, String const & parent_end_name) -> int {
		 return ptr.add_state(state, parent_index, parent_end_name); } )
		.def("add_connection",[] (MotifStateGraph  & ptr, int i, int j, String const & i_bp_name, String const & j_bp_name) {
		ptr.add_connection(i, j, i_bp_name, j_bp_name); } )
		.def("replace_state",[] (MotifStateGraph  & ptr, int i, motif::MotifStateOP const & new_state) {
		ptr.replace_state(i, new_state); } )
		.def("remove_state",[] (MotifStateGraph  & ptr, int pos) {
		ptr.remove_state(pos); } )
		.def("remove_level",[] (MotifStateGraph  & ptr, int level) {
		ptr.remove_level(level); } )
		.def("last_node",[] (MotifStateGraph  & ptr) ->  data_structure::graph::GraphNodeOP<MSNodeDataOP> {
		 return ptr.last_node(); } )
		.def("get_node",[] (MotifStateGraph const & ptr, int i) ->  data_structure::graph::GraphNodeOP<MSNodeDataOP> const & {
		 return ptr.get_node(i); } )
		.def("get_node",[] (MotifStateGraph const & ptr, util::Uuid const & uuid) ->  data_structure::graph::GraphNodeOP<MSNodeDataOP> const {
		 return ptr.get_node(uuid); } )
		.def("get_node",[] (MotifStateGraph const & ptr, String const & m_name) ->  data_structure::graph::GraphNodeOP<MSNodeDataOP> const {
		 return ptr.get_node(m_name); } )
		.def("increase_level",[] (MotifStateGraph  & ptr) {
		ptr.increase_level(); } )
		.def("to_motif_graph",[] (MotifStateGraph  & ptr) -> MotifGraphOP {
		 return ptr.to_motif_graph(); } )
		.def("unaligned_nodes",[] (MotifStateGraph const & ptr) -> data_structure::graph::GraphNodeOPs<MSNodeDataOP> const {
		 return ptr.unaligned_nodes(); } )
		.def("get_int_option",[] (MotifStateGraph  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (MotifStateGraph  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (MotifStateGraph  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (MotifStateGraph  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		 .def("set_option_value",[] (MotifGraph  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
		 .def("set_option_value",[] (MotifGraph  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
		 .def("set_option_value",[] (MotifGraph  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
		 .def("set_option_value",[] (MotifGraph  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
		;

        py::class_<MotifStateTree, std::shared_ptr<MotifStateTree>>(m, "MotifStateTree")
		// ctors
		.def(py::init<>())
		.def(py::init<MotifTreeOP const>())
		.def(py::init<MotifStateTree const>())
		.def(py::init<String const>())
		// methods
		.def("begin",[] (MotifStateTree  & ptr) -> MotifStateTree::iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateTree  & ptr) -> MotifStateTree::iterator {
		 return ptr.end(); } )
		.def("begin",[] (MotifStateTree const & ptr) -> MotifStateTree::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifStateTree const & ptr) -> MotifStateTree::const_iterator {
		 return ptr.end(); } )
		.def("add_state",[] (MotifStateTree  & ptr, motif::MotifStateOP const & state, int parent_index, int parent_end_index) -> int {
		 return ptr.add_state(state, parent_index, parent_end_index); } )
		.def("add_state",[] (MotifStateTree  & ptr, motif::MotifStateOP const & state, int parent_index, String const & parent_end_name) -> int {
		 return ptr.add_state(state, parent_index, parent_end_name); } )
		.def("add_mst",[] (MotifStateTree  & ptr, MotifStateTreeOP const & mst, int parent_index, int parent_end_index) -> int {
		 return ptr.add_mst(mst, parent_index, parent_end_index); } )
		.def("add_mst",[] (MotifStateTree  & ptr, MotifStateTreeOP const & mst, int parent_index, int parent_end_index) -> int {
		 return ptr.add_mst(mst, parent_index, parent_end_index); } )
		.def("add_connection",[] (MotifStateTree  & ptr, int i, int j, String const & i_bp_name, String const & j_bp_name) {
		ptr.add_connection(i, j, i_bp_name, j_bp_name); } )
		.def("replace_state",[] (MotifStateTree  & ptr, int i, motif::MotifStateOP const & new_state) {
		ptr.replace_state(i, new_state); } )
		.def("remove_node",[] (MotifStateTree  & ptr, int i) {
		ptr.remove_node(i); } )
		.def("remove_node_level",[] (MotifStateTree  & ptr, int level) {
		ptr.remove_node_level(level); } )
		.def("to_motif_tree",[] (MotifStateTree  & ptr) -> MotifTreeOP {
		 return ptr.to_motif_tree(); } )
		.def("topology_to_str",[] (MotifStateTree  & ptr) -> String {
		 return ptr.topology_to_str(); } )
		.def("centers",[] (MotifStateTree  & ptr) ->  math::Points {
		 return ptr.centers(); } )
		.def("connections",[] (MotifStateTree  & ptr) ->  MotifConnections const & {
		 return ptr.connections(); } )
		.def("write_pdbs",[] (MotifStateTree  & ptr, String const & fname) {
		ptr.write_pdbs(fname); } )
		.def("get_structure",[] (MotifStateTree  & ptr) ->  structure::RNAStructureOP {
		 return ptr.get_structure(); } )
		.def("size",[] (MotifStateTree  & ptr) -> size_t {
		 return ptr.size(); } )
		.def("last_node",[] (MotifStateTree  & ptr) -> MotifStateTreeNodeOP const & {
		 return ptr.last_node(); } )
		.def("get_node",[] (MotifStateTree  & ptr, int i) -> MotifStateTreeNodeOP const & {
		 return ptr.get_node(i); } )
		.def("get_node",[] (MotifStateTree  & ptr, util::Uuid const & uuid) ->  MotifStateTreeNodeOP const & {
		 return ptr.get_node(uuid); } )
		.def("get_node",[] (MotifStateTree  & ptr, String const & m_name) ->  MotifStateTreeNodeOP {
		 return ptr.get_node(m_name); } )
		.def("increase_level",[] (MotifStateTree  & ptr) {
		ptr.increase_level(); } )
		.def("decrease_level",[] (MotifStateTree  & ptr) {
		ptr.decrease_level(); } )
		.def("get_int_option",[] (MotifStateTree  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (MotifStateTree  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (MotifStateTree  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (MotifStateTree  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
        .def("set_option_value",[] (MotifGraph  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
        .def("set_option_value",[] (MotifGraph  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
        .def("set_option_value",[] (MotifGraph  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
        .def("set_option_value",[] (MotifGraph  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )

		;

        py::class_<MotifTree, std::shared_ptr<MotifTree>>(m, "MotifTree")
		// ctors
		.def(py::init<>())
		.def(py::init<String const>())
		.def(py::init<String const,MotifTreeStringType>())
		.def(py::init<MotifTree const>())
		// methods
		.def("begin",[] (MotifTree  & ptr) -> MotifTree::iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifTree  & ptr) -> MotifTree::iterator {
		 return ptr.end(); } )
		.def("begin",[] (MotifTree const & ptr) -> MotifTree::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (MotifTree const & ptr) -> MotifTree::const_iterator {
		 return ptr.end(); } )
		.def("add_motif",[] (MotifTree  & ptr, motif::MotifOP const & m, int parent_index, int parent_end_index) -> int {
		 return ptr.add_motif(m, parent_index, parent_end_index); } )
		.def("add_motif",[] (MotifTree  & ptr, motif::MotifOP const & m, int parent_index, String parent_end_name) -> int {
		 return ptr.add_motif(m, parent_index, parent_end_name); } )
		.def("add_connection",[] (MotifTree  & ptr, int i, int j, String const & i_bp_name, String const & j_bp_name) {
		ptr.add_connection(i, j, i_bp_name, j_bp_name); } )
		.def("remove_node",[] (MotifTree  & ptr, int i) {
		ptr.remove_node(i); } )
		.def("remove_node_level",[] (MotifTree  & ptr, int level) {
		ptr.remove_node_level(level); } )
		.def("to_pdb",[] (MotifTree  & ptr, String const fname, int renumber, int close_chains, int conect_statements) {
		ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
		.def("topology_to_str",[] (MotifTree  & ptr) -> String {
		 return ptr.topology_to_str(); } )
		.def("to_pretty_str",[] (MotifTree  & ptr) -> String {
		 return ptr.to_pretty_str(); } )
		.def("to_str",[] (MotifTree  & ptr) -> String {
		 return ptr.to_str(); } )
		.def("_update_merger",[] (MotifTree  & ptr) {
		ptr._update_merger(); } )
		.def("connections",[] (MotifTree  & ptr) -> MotifConnections const & {
		 return ptr.connections(); } )
		.def("beads",[] (MotifTree  & ptr) -> structure::Beads {
		 return ptr.beads(); } )
		.def("size",[] (MotifTree  & ptr) -> size_t {
		 return ptr.size(); } )
		.def("get_node",[] (MotifTree  & ptr, int i) ->  data_structure::tree::TreeNodeOP<motif::MotifOP> const & {
		 return ptr.get_node(i); } )
		.def("get_node",[] (MotifTree  & ptr, util::Uuid const & uuid) ->  data_structure::tree::TreeNodeOP<motif::MotifOP> const & {
		 return ptr.get_node(uuid); } )
		.def("get_node",[] (MotifTree  & ptr, String const & m_name) ->  data_structure::tree::TreeNodeOP<motif::MotifOP> {
		 return ptr.get_node(m_name); } )
		.def("last_node",[] (MotifTree  & ptr) ->  data_structure::tree::TreeNodeOP<motif::MotifOP> const & {
		 return ptr.last_node(); } )
		.def("write_pdbs",[] (MotifTree  & ptr, String const & fname) {
		ptr.write_pdbs(fname); } )
		.def("increase_level",[] (MotifTree  & ptr) {
		ptr.increase_level(); } )
		.def("decrease_level",[] (MotifTree  & ptr) {
		ptr.decrease_level(); } )
		.def("get_structure",[] (MotifTree  & ptr) ->  structure::RNAStructureOP const & {
		 return ptr.get_structure(); } )
		.def("secondary_structure",[] (MotifTree  & ptr) -> secondary_structure::PoseOP {
		 return ptr.secondary_structure(); } )
		.def("designable_secondary_structure",[] (MotifTree  & ptr) -> secondary_structure::PoseOP {
		 return ptr.designable_secondary_structure(); } )
		.def("get_int_option",[] (MotifTree  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (MotifTree  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (MotifTree  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (MotifTree  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
         .def("set_option_value",[] (MotifGraph  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
         .def("set_option_value",[] (MotifGraph  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
         .def("set_option_value",[] (MotifGraph  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
         .def("set_option_value",[] (MotifGraph  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
		;
/*
        py::class_<MotifTreePrinter, std::shared_ptr<MotifTreePrinter>>(m, "MotifTreePrinter")
		// methods
		.def("print_tree",[] (MotifTreePrinter  & ptr, MotifTree const & mt) -> String {
		 return ptr.print_tree(mt); } )
		;
        py::class_<_GraphtoTreeNode, std::shared_ptr<_GraphtoTreeNode>>(m, "_GraphtoTreeNode")
		// public attributes
		.def_readwrite("parent", &_GraphtoTreeNode::parent)
		.def_readwrite("parent_end_index", &_GraphtoTreeNode::parent_end_index)
		.def_readwrite("node", &_GraphtoTreeNode::node)
		.def_readwrite("motif", &_GraphtoTreeNode::motif)
		;

*/
	// enums
	py::enum_<MotifGraphStringType>(m,"MotifGraphStringType")
		.value("OLD",MotifGraphStringType::OLD)
		.value("MG",MotifGraphStringType::MG)
		.value("TOP",MotifGraphStringType::TOP)
	;
	py::enum_<MotifTreeStringType>(m,"MotifTreeStringType")
		.value("MT_STR",MotifTreeStringType::MT_STR)
		.value("TOP_STR",MotifTreeStringType::TOP_STR)
	;


}