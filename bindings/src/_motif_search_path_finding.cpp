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
using namespace motif_search::path_finding;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_motif_search_path_finding,m) {

	// free functions

    m.def("default_selector", [] () -> SelectorOP {
		 return default_selector(); } 
    );
	// classes

        py::class_<AstarScorer, std::shared_ptr<AstarScorer>>(m, "AstarScorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (AstarScorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("score",[] (AstarScorer  & ptr, Node const & node) ->  float {
		 return ptr.score(node); } )
		.def("score",[] (AstarScorer  & ptr, motif::MotifState & ms, Node const & node) ->  float {
		 return ptr.score(ms, node); } )
		// inherited methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("set_target",[] (Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
		ptr.set_target(target, target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, Node const & node) ->  float {
		 return ptr.score(node); } )
		.def("score",[] (Scorer  & ptr, motif::MotifState & state, Node const & node) ->  float {
		 return ptr.score(state, node); } )
		.def("accept_score",[] (Scorer  & ptr, Node const & node) ->  float {
		 return ptr.accept_score(node); } )
		.def("set_dummy",[] (Scorer  & ptr, float dummy) {
		ptr.set_dummy(dummy); } )
		;

        py::class_<GreedyScorer, std::shared_ptr<GreedyScorer>>(m, "GreedyScorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (GreedyScorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("score",[] (GreedyScorer  & ptr, Node const & node) ->  float {
		 return ptr.score(node); } )
		.def("score",[] (GreedyScorer  & ptr, motif::MotifState & ms, Node const & node) ->  float {
		 return ptr.score(ms, node); } )
		// inherited methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("set_target",[] (Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
		ptr.set_target(target, target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, Node const & node) ->  float {
		 return ptr.score(node); } )
		.def("score",[] (Scorer  & ptr, motif::MotifState & state, Node const & node) ->  float {
		 return ptr.score(state, node); } )
		.def("accept_score",[] (Scorer  & ptr, Node const & node) ->  float {
		 return ptr.accept_score(node); } )
		.def("set_dummy",[] (Scorer  & ptr, float dummy) {
		ptr.set_dummy(dummy); } )
		;

        py::class_<Node, std::shared_ptr<Node>>(m, "Node")
		// ctors
		.def(py::init<motif::MotifStateOP,motif_search::path_finding::NodeOP,float,int,int,int>())
		// methods
		.def("level",[] (Node const & ptr) -> int {
		 return ptr.level(); } )
		.def("size",[] (Node const & ptr) ->  int {
		 return ptr.size(); } )
		.def("node_type",[] (Node const & ptr) ->  int {
		 return ptr.node_type(); } )
		.def("ss_score",[] (Node const & ptr) ->  float {
		 return ptr.ss_score(); } )
		.def("score",[] (Node const & ptr) ->  float {
		 return ptr.score(); } )
		.def("parent_end_index",[] (Node const & ptr) ->  int {
		 return ptr.parent_end_index(); } )
		.def("state",[] (Node const & ptr) ->  motif::MotifStateOP {
		 return ptr.state(); } )
		.def("parent",[] (Node const & ptr) ->  NodeOP {
		 return ptr.parent(); } )
		;
/*
        py::class_<NodeCompare, std::shared_ptr<NodeCompare>>(m, "NodeCompare")
		// operators
		.def(py::self () motif_search::path_finding::NodeOP)
		;
*/


        py::class_<RoundRobinSelector, std::shared_ptr<RoundRobinSelector>>(m, "RoundRobinSelector")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (RoundRobinSelector const & ptr) -> Selector * {
		 return ptr.clone(); } )
		.def("add",[] (RoundRobinSelector  & ptr, String const & lib_name) {
		ptr.add(lib_name); } )
		.def("add",[] (RoundRobinSelector  & ptr, String const & lib_name, motif::MotifStateEnsembleOP mse) {
		ptr.add(lib_name, mse); } )
		.def("add",[] (RoundRobinSelector  & ptr, motif::MotifOP motif) {
		ptr.add(motif); } )
		// inherited methods
		.def("clone",[] (Selector const & ptr) -> Selector * {
		 return ptr.clone(); } )
		.def("add",[] (Selector  & ptr, String const & lib_name) {
		ptr.add(lib_name); } )
		.def("add",[] (Selector  & ptr, String const & lib_name, motif::MotifStateEnsembleOP mse) {
		ptr.add(lib_name, mse); } )
		.def("add",[] (Selector  & ptr, motif::MotifOP motif) {
		ptr.add(motif); } )
		.def("connect",[] (Selector  & ptr, String const & name_i, String const & name_j) {
		ptr.connect(name_i, name_j); } )
		.def("start",[] (Selector const & ptr, int parent_type) {
		ptr.start(parent_type); } )
		.def("next",[] (Selector const & ptr) -> SelectorNodeDataOP {
		 return ptr.next(); } )
		.def("finished",[] (Selector const & ptr) -> bool {
		 return ptr.finished(); } )
		.def("size",[] (Selector const & ptr) -> size_t {
		 return ptr.size(); } )
		;
/*
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("set_target",[] (Scorer  & ptr, motif_search::path_findingstructure::BasepairStateOP target, bool target_an_aligned_end) {
		ptr.set_target(target, target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, Node const &) ->  float {
		 return ptr.score(&); } )
		.def("score",[] (Scorer  & ptr, motif::MotifState &, Node const &) ->  float {
		 return ptr.score(&, &); } )
		.def("accept_score",[] (Scorer  & ptr, Node const & node) ->  float {
		 return ptr.accept_score(node); } )
		.def("set_dummy",[] (Scorer  & ptr, float dummy) {
		ptr.set_dummy(dummy); } )
		;
*/
        py::class_<Search, std::shared_ptr<Search>>(m, "Search")
		// ctors
		.def(py::init<ScorerOP,SelectorOP,motif_search::SolutionFilterOP>())
		// methods
		.def("clone",[] (Search const & ptr) -> motif_search::Search * {
		 return ptr.clone(); } )
		.def("setup",[] (Search  & ptr, motif_search:: ProblemOP prblm) {
		ptr.setup(prblm); } )
		.def("start",[] (Search  & ptr) {
		ptr.start(); } )
		.def("finished",[] (Search  & ptr) -> bool {
		 return ptr.finished(); } )
		.def("next",[] (Search  & ptr) -> motif_search::SolutionOP {
		 return ptr.next(); } )
		;

        py::class_<Selector, std::shared_ptr<Selector>>(m, "Selector")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Selector const & ptr) -> Selector * {
		 return ptr.clone(); } )
		.def("add",[] (Selector  & ptr, String const & lib_name) {
		ptr.add(lib_name); } )
		.def("add",[] (Selector  & ptr, String const & lib_name, motif::MotifStateEnsembleOP mse) {
		ptr.add(lib_name, mse); } )
		.def("add",[] (Selector  & ptr, motif::MotifOP motif) {
		ptr.add(motif); } )
		.def("connect",[] (Selector  & ptr, String const & name_i, String const & name_j) {
		ptr.connect(name_i, name_j); } )
		.def("start",[] (Selector const & ptr, int parent_type) {
		ptr.start(parent_type); } )
		.def("next",[] (Selector const & ptr) -> SelectorNodeDataOP {
		 return ptr.next(); } )
		.def("finished",[] (Selector const & ptr) -> bool {
		 return ptr.finished(); } )
		.def("size",[] (Selector const & ptr) -> size_t {
		 return ptr.size(); } )
		;

        py::class_<SelectorNodeData, std::shared_ptr<SelectorNodeData>>(m, "SelectorNodeData")
		// ctors
		.def(py::init<String const &,motif::MotifStateOPs const &,int>())
		// public attributes
		.def_readwrite("name", &SelectorNodeData::name)
		.def_readwrite("motif_states", &SelectorNodeData::motif_states)
		.def_readwrite("type", &SelectorNodeData::type)
		;



}