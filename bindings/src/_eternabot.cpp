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
using namespace eternabot;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_eternabot,m) {

	// classes

        py::class_<ABasicTest, std::shared_ptr<ABasicTest>>(m, "ABasicTest")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (ABasicTest  & ptr, FeaturesOP const & features) ->  float {
		 return ptr.score(features); } )
		// inherited methods

		.def("mean",[] (Strategy const & ptr) ->  float {
		 return ptr.mean(); } )
		.def("stdev",[] (Strategy const & ptr) ->  float {
		 return ptr.stdev(); } )
		;

        py::class_<BerexTest, std::shared_ptr<BerexTest>>(m, "BerexTest")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (BerexTest  & ptr, FeaturesOP const & features) -> float {
		 return ptr.score(features); } )
		// inherited methods
		.def("mean",[] (Strategy const & ptr) ->  float {
		 return ptr.mean(); } )
		.def("stdev",[] (Strategy const & ptr) ->  float {
		 return ptr.stdev(); } )
		;

        py::class_<CleanPlotStackCapsandSafeGC, std::shared_ptr<CleanPlotStackCapsandSafeGC>>(m, "CleanPlotStackCapsandSafeGC")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (CleanPlotStackCapsandSafeGC  & ptr, FeaturesOP const & features) -> float {
		 return ptr.score(features); } )
		// inherited methods
		.def("mean",[] (Strategy const & ptr) ->  float {
		 return ptr.mean(); } )
		.def("stdev",[] (Strategy const & ptr) ->  float {
		 return ptr.stdev(); } )
		;

        py::class_<DirectionofGCPairsinMultiLoops, std::shared_ptr<DirectionofGCPairsinMultiLoops>>(m, "DirectionofGCPairsinMultiLoops")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (DirectionofGCPairsinMultiLoops  & ptr, FeaturesOP const & features) -> float {
		 return ptr.score(features); } )
		// inherited methods
		.def("mean",[] (Strategy const & ptr) ->  float {
		 return ptr.mean(); } )
		.def("stdev",[] (Strategy const & ptr) ->  float {
		 return ptr.stdev(); } )
		;

        py::class_<FeatureGenerator, std::shared_ptr<FeatureGenerator>>(m, "FeatureGenerator")
		// ctors
		.def(py::init<>())
		// methods
		.def("get_features",[] (FeatureGenerator  & ptr, secondary_structure::PoseOP const & p) -> FeaturesOP {
		 return ptr.get_features(p); } )
		.def("update_features",[] (FeatureGenerator  & ptr, FeaturesOP & features, secondary_structure::PoseOP const & p) {
		ptr.update_features(features, p); } )
		;

        py::class_<Features, std::shared_ptr<Features>>(m, "Features")
		// ctors
		.def(py::init<>())
		// public attributes
		.def_readwrite("length", &Features::length)
		.def_readwrite("a_count", &Features::a_count)
		.def_readwrite("c_count", &Features::c_count)
		.def_readwrite("g_count", &Features::g_count)
		.def_readwrite("u_count", &Features::u_count)
		.def_readwrite("gu", &Features::gu)
		.def_readwrite("gc", &Features::gc)
		.def_readwrite("ua", &Features::ua)
		.def_readwrite("meltpoint", &Features::meltpoint)
		.def_readwrite("fe", &Features::fe)
		.def_readwrite("dotplot", &Features::dotplot)
		.def_readwrite("pairmap", &Features::pairmap)
		.def_readwrite("helices", &Features::helices)
		;

        py::class_<NumofYellowNucleotidesperLengthofString, std::shared_ptr<NumofYellowNucleotidesperLengthofString>>(m, "NumofYellowNucleotidesperLengthofString")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (NumofYellowNucleotidesperLengthofString  & ptr, FeaturesOP const & features) -> float {
		 return ptr.score(features); } )
		// inherited methods
		.def("mean",[] (Strategy const & ptr) ->  float {
		 return ptr.mean(); } )
		.def("stdev",[] (Strategy const & ptr) ->  float {
		 return ptr.stdev(); } )
		;
/*
        py::class_<Parameters, std::shared_ptr<Parameters>>(m, "Parameters")
		;
*/
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("setup",[] (Scorer  & ptr, secondary_structure::PoseOP const & p) {
		ptr.setup(p); } )
		.def("score_secondary_structure",[] (Scorer  & ptr, secondary_structure::PoseOP const & p) -> float {
		 return ptr.score_secondary_structure(p); } )
		.def("scores",[] (Scorer  & ptr) -> Floats const & {
		 return ptr.scores(); } )
		;

        py::class_<SequenceDesigner, std::shared_ptr<SequenceDesigner>>(m, "SequenceDesigner")
		// ctors
		.def(py::init<>())
		// methods
		.def("setup",[] (SequenceDesigner  & ptr) {
		ptr.setup(); } )
		.def("design",[] (SequenceDesigner  & ptr, secondary_structure::PoseOP const & p) -> SequenceDesignerResultOPs const &  {
		 return ptr.design(p); } )
		.def("get_int_option",[] (SequenceDesigner  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (SequenceDesigner  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (SequenceDesigner  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (SequenceDesigner  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		.def("has_option",[] (SequenceDesigner  & ptr, String const & name) ->  bool {
		 return ptr.has_option(name); } )
		.def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
		;

        py::class_<SequenceDesignerResult, std::shared_ptr<SequenceDesignerResult>>(m, "SequenceDesignerResult")
		// ctors
		.def(py::init<String const &,float>())
		// public attributes
		.def_readwrite("sequence", &SequenceDesignerResult::sequence)
		.def_readwrite("score", &SequenceDesignerResult::score)
		;
/*
        py::class_<Strategy, std::shared_ptr<Strategy>>(m, "Strategy")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (Strategy  & ptr, FeaturesOP const &) -> float {
		 return ptr.score(&); } )
		.def("mean",[] (Strategy const & ptr) ->  float {
		 return ptr.mean(); } )
		.def("stdev",[] (Strategy const & ptr) ->  float {
		 return ptr.stdev(); } )
		;
        py::class_<sequence_designer_result_less_than_key, std::shared_ptr<sequence_designer_result_less_than_key>>(m, "sequence_designer_result_less_than_key")
		// operators
		.def(py::self () eternabot::SequenceDesignerResultOP)
		;

*/


}