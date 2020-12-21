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
using namespace sequence_optimization;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_sequence_optimization,m) {

	// classes
/*
        py::class_<DesignableBP, std::shared_ptr<DesignableBP>>(m, "DesignableBP")
		// methods
		.def("m_id_bot",[] (DesignableBP  & ptr, nullptr) , m_id_top(nullptr ) -> ), {
		 return ptr.m_id_bot(, ); } )
		.def("update_state",[] (DesignableBP  & ptr, Strings const & bp_name) {
		ptr.update_state(bp_name); } )
		.def("revert_state",[] (DesignableBP  & ptr) {
		ptr.revert_state(); } )
		// public attributes
		.def_readwrite("bp", &DesignableBP::bp)
		.def_readwrite("last_state", &DesignableBP::last_state)
		.def_readwrite("m_id_bot", &DesignableBP::m_id_bot)
		.def_readwrite("m_id_top", &DesignableBP::m_id_top)
		;
*/
        py::class_<ExternalTargetScorer, std::shared_ptr<ExternalTargetScorer>>(m, "ExternalTargetScorer")
		// ctors
		.def(py::init<structure::BasepairStateOP const &,int,int,bool>())
		// methods
		.def("score",[] (ExternalTargetScorer  & ptr, motif_data_structure::MotifStateGraphOP const & msg) -> float {
		 return ptr.score(msg); } )
		// inherited methods
		.def("score",[] (SequenceOptimizerScorer  & ptr, motif_data_structure::MotifStateGraphOP const & graph) -> float {
		 return ptr.score(graph); } )
		.def("motif_state_diff",[] (SequenceOptimizerScorer  & ptr, structure::BasepairStateOP const & end1, structure::BasepairStateOP const & end2) -> float {
		 return ptr.motif_state_diff(end1, end2); } )
		;

        py::class_<InternalTargetScorer, std::shared_ptr<InternalTargetScorer>>(m, "InternalTargetScorer")
		// ctors
		.def(py::init<int,int,int,int,bool>())
		// methods
		.def("score",[] (InternalTargetScorer  & ptr, motif_data_structure::MotifStateGraphOP const & msg) -> float {
		 return ptr.score(msg); } )
		// inherited methods
		.def("score",[] (SequenceOptimizerScorer  & ptr, motif_data_structure::MotifStateGraphOP const & graph) -> float {
		 return ptr.score(graph); } )
		.def("motif_state_diff",[] (SequenceOptimizerScorer  & ptr, structure::BasepairStateOP const & end1, structure::BasepairStateOP const & end2) -> float {
		 return ptr.motif_state_diff(end1, end2); } )
		;

        py::class_<MultiTargetScorer, std::shared_ptr<MultiTargetScorer>>(m, "MultiTargetScorer")
		// ctors
		.def(py::init<std::vector<SequenceOptimizerScorerOP>>())
		// methods
		.def("score",[] (MultiTargetScorer  & ptr, motif_data_structure::MotifStateGraphOP const & msg) -> float {
		 return ptr.score(msg); } )
		// inherited methods
		.def("score",[] (SequenceOptimizerScorer  & ptr, motif_data_structure::MotifStateGraphOP const & graph) -> float {
		 return ptr.score(graph); } )
		.def("motif_state_diff",[] (SequenceOptimizerScorer  & ptr, structure::BasepairStateOP const & end1, structure::BasepairStateOP const & end2) -> float {
		 return ptr.motif_state_diff(end1, end2); } )
		;

        py::class_<OptimizedSequence, std::shared_ptr<OptimizedSequence>>(m, "OptimizedSequence")
		// public attributes
		.def_readwrite("sequence", &OptimizedSequence::sequence)
		.def_readwrite("dist_score", &OptimizedSequence::dist_score)
		.def_readwrite("eterna_score", &OptimizedSequence::eterna_score)
		;

        py::class_<SequenceOptimizer3D, std::shared_ptr<SequenceOptimizer3D>>(m, "SequenceOptimizer3D")
		// ctors
		.def(py::init<>())
		// methods
		.def("set_scorer",[] (SequenceOptimizer3D  & ptr, sequence_optimization::SequenceOptimizerScorerOP scorer) {
		ptr.set_scorer(scorer); } )
		.def("get_optimized_sequences",[] (SequenceOptimizer3D  & ptr, motif_data_structure::MotifGraphOP const & mg, SequenceOptimizerScorerOP scorer) ->  OptimizedSequenceOPs {
		 return ptr.get_optimized_sequences(mg, scorer); } )
		.def("get_optimized_sequences",[] (SequenceOptimizer3D  & ptr, motif_data_structure::MotifGraphOP const & mg) -> OptimizedSequenceOPs {
		 return ptr.get_optimized_sequences(mg); } )
		.def("get_optimized_mg",[] (SequenceOptimizer3D  & ptr, motif_data_structure::MotifGraphOP const & mg, ::SequenceOptimizerScorerOP scorer) ->  motif_data_structure::MotifGraphOP {
		 return ptr.get_optimized_mg(mg, scorer); } )
		.def("get_optimized_mg",[] (SequenceOptimizer3D  & ptr, motif_data_structure::MotifGraphOP const & mg) -> motif_data_structure::MotifGraphOP  {
		 return ptr.get_optimized_mg(mg); } )
		.def("options",[] (SequenceOptimizer3D  & ptr) ->  base::Options & {
		 return ptr.options(); } )
		.def("get_int_option",[] (SequenceOptimizer3D  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (SequenceOptimizer3D  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (SequenceOptimizer3D  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (SequenceOptimizer3D  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		.def("set_option_value",[] (SequenceOptimizer3D  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (SequenceOptimizer3D  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (SequenceOptimizer3D  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (SequenceOptimizer3D  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
		.def("update_var_options",[] (SequenceOptimizer3D  & ptr) {
		ptr.update_var_options(); } )
		;
/*
        py::class_<SequenceOptimizerScorer, std::shared_ptr<SequenceOptimizerScorer>>(m, "SequenceOptimizerScorer")
		// ctors
		.def(py::init<bool>())
		// methods
		.def("score",[] (SequenceOptimizerScorer  & ptr, motif_data_structure::MotifStateGraphOP const & graph) -> float {
		 return ptr.score(&); } )
		.def("motif_state_diff",[] (SequenceOptimizerScorer  & ptr, sequence_optimizationstructure::BasepairStateOP const & end1, sequence_optimizationstructure::BasepairStateOP const & end2) -> float {
		 return ptr.motif_state_diff(end1, end2); } )
		;
*/


}