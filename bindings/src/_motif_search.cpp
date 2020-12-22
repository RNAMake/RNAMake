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
using namespace motif_search;
using namespace motif_search::path_finding;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_motif_search,m) {

	// classes

        py::class_<MotifStateMonteCarlo, std::shared_ptr<MotifStateMonteCarlo>>(m, "MotifStateMonteCarlo")
		// ctors
		.def(py::init<std::vector<motif::MotifStateOPs> const &>())
		// methods
		.def("setup",[] (MotifStateMonteCarlo  & ptr, motif_data_structure::MotifStateGraphOP msg, int arg1, int arg2, int arg3, int arg4, bool arg5) {
		ptr.setup(msg, arg1, arg2, arg3, arg4, arg5); } )
		.def("run",[] (MotifStateMonteCarlo  & ptr) {
		ptr.run(); } )
		.def("start",[] (MotifStateMonteCarlo  & ptr) {
		ptr.start(); } )
		.def("next",[] (MotifStateMonteCarlo  & ptr) -> MotifStateMonteCarloSolutionOP {
		 return ptr.next(); } )
		.def("next_state",[] (MotifStateMonteCarlo  & ptr) -> MotifStateMonteCarloSolutionNewOP {
		 return ptr.next_state(); } )
		.def("finished",[] (MotifStateMonteCarlo  & ptr) -> bool {
		 return ptr.finished(); } )
		.def("lookup",[] (MotifStateMonteCarlo  & ptr, util::StericLookupNew const & sl) {
		ptr.lookup(sl); } )
		.def("options",[] (MotifStateMonteCarlo  & ptr) ->  base::Options & {
		 return ptr.options(); } )
		.def("get_int_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		.def("has_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  bool {
		 return ptr.has_option(name); } )
		.def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
		;

        py::class_<MotifStateMonteCarloSolution, std::shared_ptr<MotifStateMonteCarloSolution>>(m, "MotifStateMonteCarloSolution")
		// ctors
		.def(py::init<motif_data_structure::MotifGraphOP,float>())
		// public attributes
		.def_readwrite("mg", &MotifStateMonteCarloSolution::mg)
		.def_readwrite("score", &MotifStateMonteCarloSolution::score)
		;

        py::class_<MotifStateMonteCarloSolutionNew, std::shared_ptr<MotifStateMonteCarloSolutionNew>>(m, "MotifStateMonteCarloSolutionNew")
		// ctors
		.def(py::init<motif_data_structure::MotifStateGraphOP,float>())
		// public attributes
		.def_readwrite("msg", &MotifStateMonteCarloSolutionNew::msg)
		.def_readwrite("score", &MotifStateMonteCarloSolutionNew::score)
		;

        py::class_<NoExclusionFilter, std::shared_ptr<NoExclusionFilter>>(m, "NoExclusionFilter")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (NoExclusionFilter const & ptr) -> SolutionFilter * {
		 return ptr.clone(); } )
		.def("accept",[] (NoExclusionFilter  & ptr, Strings const & motif_names) -> bool {
		 return ptr.accept(motif_names); } )
		// inherited methods
		.def("clone",[] (SolutionFilter const & ptr) -> SolutionFilter * {
		 return ptr.clone(); } )
		.def("accept",[] (SolutionFilter  & ptr, Strings const & sols) -> bool {
		 return ptr.accept(sols); } )
		;

/*
        py::class_<Node, std::shared_ptr<Node>>(m, "Node")
		// methods
		.def("get_type",[] (Node const & ptr) ->  NodeType {
		 return ptr.get_type(); } )
		.def("get_lib_name",[] (Node const & ptr) ->  String const & {
		 return ptr.get_lib_name(); } )
		.def("get_motif_state",[] (Node const & ptr) ->  motif::MotifStateOP {
		 return ptr.get_motif_state(); } )
		.def("get_motif_state_ensemble",[] (Node const & ptr) ->  motif::MotifStateEnsembleOP {
		 return ptr.get_motif_state_ensemble(); } )
		;
        py::enum_<NodeType, std::shared_ptr<NodeType>>(m, "NodeType")

		;
        py::class_<path_finding::Parameters, std::shared_ptr<path_finding::Parameters>>(m, "Parameters")
		// public attributes
		.def_readwrite("max_helix_size", &Parameters::max_helix_size)
		.def_readwrite("min_helix_size", &path_finding::Parameters::min_helix_size)
		;
*/
    py::class_<motif_search::path_finding::Parameters, std::shared_ptr<motif_search::path_finding::Parameters>>(m, "motif_search::path_finding::Parameters")
            // public attributes
            .def_readwrite("sterics", &motif_search::path_finding::Parameters::sterics)
            .def_readwrite("helix_end", &motif_search::path_finding::Parameters::helix_end)
            .def_readwrite("max_node_level", &motif_search::path_finding::Parameters::max_node_level)
            .def_readwrite("min_size", &motif_search::path_finding::Parameters::min_size)
            .def_readwrite("max_size", &motif_search::path_finding::Parameters::max_size)
            .def_readwrite("max_solutions", &motif_search::path_finding::Parameters::max_solutions)
            .def_readwrite("min_node_level", &motif_search::path_finding::Parameters::min_node_level)
            .def_readwrite("accept_score", &motif_search::path_finding::Parameters::accept_score)
            .def_readwrite("min_ss_score", &motif_search::path_finding::Parameters::min_ss_score)
            .def_readwrite("max_steps", &motif_search::path_finding::Parameters::max_steps)
            .def_readwrite("return_best", &motif_search::path_finding::Parameters::return_best)
            ;
        py::class_<Problem, std::shared_ptr<Problem>>(m, "Problem")
		// ctors
		.def(py::init<structure::BasepairStateOP,structure::BasepairStateOP,util::StericLookupNewOP,bool>())
		// public attributes
		.def_readwrite("start", &Problem::start)
		.def_readwrite("end", &Problem::end)
		.def_readwrite("lookup", &Problem::lookup)
		.def_readwrite("target_an_aligned_end", &Problem::target_an_aligned_end)
		;

        py::class_<RemoveDuplicateHelices, std::shared_ptr<RemoveDuplicateHelices>>(m, "RemoveDuplicateHelices")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (RemoveDuplicateHelices const & ptr) -> SolutionFilter * {
		 return ptr.clone(); } )
		.def("accept",[] (RemoveDuplicateHelices  & ptr, Strings const & motif_names) -> bool {
		 return ptr.accept(motif_names); } )
		// inherited methods
		.def("clone",[] (SolutionFilter const & ptr) -> SolutionFilter * {
		 return ptr.clone(); } )
		.def("accept",[] (SolutionFilter  & ptr, Strings const & strs) -> bool {
		 return ptr.accept(strs); } )
		;
/*
        py::class_<motif_search::Search, std::shared_ptr<motif_search::Search>>(m, "motif_search::Search")
		// ctors
		.def(py::init<String const &>())
		// methods
		.def("clone",[] (motif_search::Search const & ptr) -> motif_search::Search * {
		 return ptr.clone(); } )
		.def("setup",[] (motif_search::Search  & ptr, motif_search::ProblemOP ) {
		ptr.setup(); } )
		.def("start",[] (motif_search::Search  & ptr) {
		ptr.start(); } )
		.def("finished",[] (motif_search::Search  & ptr) -> bool {
		 return ptr.finished(); } )
		.def("next",[] (motif_search::Search  & ptr) -> SolutionOP {
		 return ptr.next(); } )
		.def("name",[] (motif_search::Search  & ptr) -> String const & {
		 return ptr.name(); } )
		.def("get_int_option",[] (motif_search::Search  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (motif_search::Search  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (motif_search::Search  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (motif_search::Search  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		.def("set_option_value",[] (motif_search::Search  & ptr, String const & name, T const & val) {
		ptr.set_option_value(name, val); } )
		;
*/
        py::class_<Solution, std::shared_ptr<Solution>>(m, "Solution")
		// ctors
		.def(py::init<motif_data_structure::MotifStateGraphOP,float>())
		// public attributes
		.def_readwrite("graph", &Solution::graph)
		.def_readwrite("score", &Solution::score)
		;
/*
        py::class_<SolutionFilter, std::shared_ptr<SolutionFilter>>(m, "SolutionFilter")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (SolutionFilter const & ptr) -> SolutionFilter * {
		 return ptr.clone(); } )
		.def("accept",[] (SolutionFilter  & ptr, Strings const &) -> bool {
		 return ptr.accept(&); } )
		;
*/
        py::class_<SolutionToplogy, std::shared_ptr<SolutionToplogy>>(m, "SolutionToplogy")
		// ctors
		.def(py::init<motif_data_structure::MotifStateEnsembleOPGraphOP>())
		// methods
		.def("begin",[] (SolutionToplogy  & ptr) -> SolutionToplogy::iterator {
		 return ptr.begin(); } )
		.def("end",[] (SolutionToplogy  & ptr) -> SolutionToplogy::iterator {
		 return ptr.end(); } )
		.def("begin",[] (SolutionToplogy const & ptr) -> SolutionToplogy::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (SolutionToplogy const & ptr) -> SolutionToplogy::const_iterator {
		 return ptr.end(); } )
		.def("initialize_solution",[] (SolutionToplogy  & ptr, structure::BasepairStateOP bp_state) -> motif_data_structure::MotifStateGraphOP {
		 return ptr.initialize_solution(bp_state); } )
		.def("initialize_solution_no_start",[] (SolutionToplogy  & ptr, structure::BasepairStateOP bp_state) -> motif_data_structure::MotifStateGraphOP {
		 return ptr.initialize_solution_no_start(bp_state); } )
		.def("get_motif_state",[] (SolutionToplogy  & ptr, Index pos) -> motif::MotifStateOP {
		 return ptr.get_motif_state(pos); } )
		.def("get_solution_nie",[] (SolutionToplogy  & ptr) -> std::vector<data_structure::NodeIndexandEdge> const & {
		 return ptr.get_solution_nie(); } )
		.def("size",[] (SolutionToplogy  & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("get_ensemble_size",[] (SolutionToplogy  & ptr, Index pos) ->  size_t {
		 return ptr.get_ensemble_size(pos); } )
		;

        py::class_<SolutionToplogyFactory, std::shared_ptr<SolutionToplogyFactory>>(m, "SolutionToplogyFactory")
		// ctors
		.def(py::init<>())
		// methods
		.def("generate_toplogy",[] (SolutionToplogyFactory  & ptr, SolutionTopologyTemplate const & sol_template) -> SolutionToplogyOP {
		 return ptr.generate_toplogy(sol_template); } )
		.def("get_int_option",[] (SolutionToplogyFactory  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (SolutionToplogyFactory  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (SolutionToplogyFactory  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (SolutionToplogyFactory  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
         .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
         .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
         .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
         .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
		;

        py::class_<SolutionTopologyTemplate, std::shared_ptr<SolutionTopologyTemplate>>(m, "SolutionTopologyTemplate")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (SolutionTopologyTemplate  & ptr) -> SolutionTopologyTemplate::iterator {
		 return ptr.begin(); } )
		.def("end",[] (SolutionTopologyTemplate  & ptr) -> SolutionTopologyTemplate::iterator {
		 return ptr.end(); } )
		.def("begin",[] (SolutionTopologyTemplate const & ptr) -> SolutionTopologyTemplate::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (SolutionTopologyTemplate const & ptr) -> SolutionTopologyTemplate::const_iterator {
		 return ptr.end(); } )
		.def("add_library",[] (SolutionTopologyTemplate  & ptr, String const & lib_name) {
		ptr.add_library(lib_name); } )
		.def("add_library",[] (SolutionTopologyTemplate  & ptr, String const & lib_name, data_structure::NodeIndexandEdge const & parent_nie) {
		ptr.add_library(lib_name, parent_nie); } )
		.def("add_motif_state",[] (SolutionTopologyTemplate  & ptr, motif::MotifStateOP ms) {
		ptr.add_motif_state(ms); } )
		.def("add_motif_state",[] (SolutionTopologyTemplate  & ptr, motif::MotifStateOP ms, data_structure::NodeIndexandEdge const & parent_nie) {
		ptr.add_motif_state(ms, parent_nie); } )
		.def("add_ensemble",[] (SolutionTopologyTemplate  & ptr, motif::MotifStateEnsembleOP mse) {
		ptr.add_ensemble(mse); } )
		.def("add_ensemble",[] (SolutionTopologyTemplate  & ptr, motif::MotifStateEnsembleOP mse, data_structure::NodeIndexandEdge const & parent_nie) {
		ptr.add_ensemble(mse, parent_nie); } )
		.def("has_parent",[] (SolutionTopologyTemplate const & ptr, Index ni) ->  bool {
		 return ptr.has_parent(ni); } )
		.def("get_parent_index",[] (SolutionTopologyTemplate const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_index(ni); } )
		.def("get_parent_end_index",[] (SolutionTopologyTemplate const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_end_index(ni); } )
		;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_search::exhaustive
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes

    py::class_<motif_search::exhaustive::DefaultScorer, std::shared_ptr<motif_search::exhaustive::DefaultScorer>>(m, "ExhaustiveDefaultScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (motif_search::exhaustive::DefaultScorer const & ptr) -> motif_search::exhaustive::Scorer * {
                return ptr.clone(); } )
            .def("score",[] (motif_search::exhaustive::DefaultScorer  & ptr, structure::BasepairState const & bps) ->  float {
                return ptr.score(bps); } )
                    // inherited methods
            .def("clone",[] (motif_search::exhaustive::Scorer const & ptr) -> motif_search::exhaustive::Scorer * {
                return ptr.clone(); } )
            .def("set_target",[] (motif_search::exhaustive::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                ptr.set_target(target, target_an_aligned_end); } )
            .def("score",[] (motif_search::exhaustive::Scorer  & ptr, structure::BasepairState const & bps) ->  float {
                return ptr.score(bps); } )
            ;

    py::class_<motif_search::exhaustive::MotifStateEnumerator, std::shared_ptr<motif_search::exhaustive::MotifStateEnumerator>>(m, "MotifStateEnumerator")
            // ctors
            .def(py::init<motif_search::SolutionToplogy>())
                    // methods
            .def("start",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr, structure:: BasepairStateOP bp) {
                ptr.start(bp); } )
            .def("finished",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr) -> bool {
                return ptr.finished(); } )
            .def("next",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr) {
                ptr.next(); } )
            .def("top_state",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr) -> motif::MotifStateOP {
                return ptr.top_state(); } )
            .def("all_states",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr) -> motif::MotifStateOPs const & {
                return ptr.all_states(); } )
            .def("set_size_limit",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr, int size_limit) {
                ptr.set_size_limit(size_limit); } )
            ;

/*
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("set_target",[] (Scorer  & ptr, motif_search::exhaustivestructure::BasepairStateOP target, bool target_an_aligned_end) {
		ptr.set_target(target, target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, structure::BasepairState const & bps) ->  float {
		 return ptr.score(bps); } )
		;
*/
    py::class_<motif_search::exhaustive::ScorerFactory, std::shared_ptr<motif_search::exhaustive::ScorerFactory>>(m, "ExhaustiveScorerFactory")
            // ctors
            .def(py::init<>())
                    // methods
            .def("get_scorer",[] (motif_search::exhaustive::ScorerFactory  & ptr, String const & scorer_name) -> motif_search::exhaustive::ScorerOP {
                return ptr.get_scorer(scorer_name); } )
            ;

    py::class_<motif_search::exhaustive::Search, std::shared_ptr<motif_search::exhaustive::Search>>(m, "ExhaustiveSearch")
            // ctors
            .def(py::init<motif_search::exhaustive::ScorerOP, motif_search::SolutionToplogy const & , motif_search::SolutionFilterOP>())
                    // methods
            .def("clone",[] (motif_search::exhaustive::Search const & ptr) -> motif_search::Search * {
                return ptr.clone(); } )
            .def("setup",[] (motif_search::exhaustive::Search  & ptr, motif_search::ProblemOP p) {
                ptr.setup(p); } )
            .def("start",[] (motif_search::exhaustive::Search  & ptr) {
                ptr.start(); } )
            .def("finished",[] (motif_search::exhaustive::Search  & ptr) -> bool {
                return ptr.finished(); } )
            .def("next",[] (motif_search::exhaustive::Search  & ptr) -> motif_search::SolutionOP {
                return ptr.next(); } )
            ;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_search::monte_carlo
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // classes

    py::class_<motif_search::monte_carlo::DefaultScorer, std::shared_ptr<motif_search::monte_carlo::DefaultScorer>>(m, "MonteCarloDefaultScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (motif_search::monte_carlo::DefaultScorer const & ptr) -> motif_search::monte_carlo::Scorer * {
                return ptr.clone(); } )
            .def("score",[] (motif_search::monte_carlo::DefaultScorer  & ptr, structure::BasepairState const & bps) ->  float {
                return ptr.score(bps); } )
                    // inherited methods
            .def("clone",[] (Scorer const & ptr) -> Scorer * {
                return ptr.clone(); } )
            .def("set_target",[] (motif_search::monte_carlo::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                ptr.set_target(target, target_an_aligned_end); } )
            .def("score",[] (motif_search::monte_carlo::Scorer  & ptr, structure::BasepairState const & bps) -> float {
                return ptr.score(bps); } )
            ;

    py::class_<motif_search::monte_carlo::MotifSwapMove, std::shared_ptr<motif_search::monte_carlo::MotifSwapMove>>(m, "MotifSwapMove")
            // ctors
            .def(py::init<motif_search::monte_carlo::ScorerOP, motif_search::SolutionToplogy const &>())
                    // methods
            .def("apply",[] (motif_search::monte_carlo::MotifSwapMove  & ptr, motif_data_structure::MotifStateGraphOP msg, float current_score) -> bool {
                return ptr.apply(msg, current_score); } )
            .def("score",[] (motif_search::monte_carlo::MotifSwapMove  & ptr) -> float {
                return ptr.score(); } )
            .def("undo",[] (motif_search::monte_carlo::MotifSwapMove  & ptr, motif_data_structure::MotifStateGraphOP msg) {
                ptr.undo(msg); } )
                    // inherited methods
            .def("apply",[] (motif_search::monte_carlo::Move  & ptr, motif_data_structure::MotifStateGraphOP msgop, float num) -> bool {
                return ptr.apply(msgop, num); } )
            .def("score",[] (motif_search::monte_carlo::Move  & ptr) -> float {
                return ptr.score(); } )
            .def("undo",[] (motif_search::monte_carlo::Move  & ptr, motif_data_structure::MotifStateGraphOP msgop) {
                ptr.undo(msgop); } )
            .def("set_temperature",[] (motif_search::monte_carlo::Move  & ptr, float temp) {
                ptr.set_temperature(temp); } )
            .def("scale_temperature",[] (motif_search::monte_carlo::Move  & ptr, float scale) {
                ptr.scale_temperature(scale); } )
            ;
/*
        py::class_<Move, std::shared_ptr<Move>>(m, "Move")
		// ctors
		.def(py::init<String const &>())
		// methods
		.def("apply",[] (Move  & ptr, motif_data_structure:: MotifStateGraphOP, float ) -> bool {
		 return ptr.apply(MotifStateGraphOP, ); } )
		.def("score",[] (Move  & ptr) -> float {
		 return ptr.score(); } )
		.def("undo",[] (Move  & ptr, motif_data_structure:: MotifStateGraphOP) {
		ptr.undo(MotifStateGraphOP); } )
		.def("set_temperature",[] (Move  & ptr, float temp) {
		ptr.set_temperature(temp); } )
		.def("scale_temperature",[] (Move  & ptr, float scale) {
		ptr.scale_temperature(scale); } )
		;
		;

        py::class_<Parameters, std::shared_ptr<Parameters>>(m, "Parameters")
		// public attributes
		.def_readwrite("accept_score", &Parameters::accept_score)
		.def_readwrite("max_size", &Parameters::max_size)
		;
*/
    py::class_< motif_search::monte_carlo::MoveSet, std::shared_ptr< motif_search::monte_carlo::MoveSet>>(m, "MoveSet");

    py::class_< motif_search::monte_carlo::ScaledScorer, std::shared_ptr< motif_search::monte_carlo::ScaledScorer>>(m, "ScaledScorer")
            // ctors
            .def(py::init<float,float>())
                    // methods
            .def("clone",[] ( motif_search::monte_carlo::ScaledScorer const & ptr) ->  motif_search::monte_carlo::Scorer * {
                return ptr.clone(); } )
            .def("score",[] ( motif_search::monte_carlo::ScaledScorer  & ptr, structure::BasepairState const & bps) ->  float {
                return ptr.score(bps); } )
                    // inherited methods
            .def("clone",[] ( motif_search::monte_carlo::Scorer const & ptr) ->  motif_search::monte_carlo::Scorer * {
                return ptr.clone(); } )
            .def("set_target",[] ( motif_search::monte_carlo::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                ptr.set_target(target, target_an_aligned_end); } )
            .def("score",[] ( motif_search::monte_carlo::Scorer  & ptr, structure::BasepairState const & bps) -> float {
                return ptr.score(bps); } )
            ;
/*
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("set_target",[] (Scorer  & ptr, motif_search::monte_carlostructure::BasepairStateOP target, bool target_an_aligned_end) {
		ptr.set_target(target, target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, structure::BasepairState const & bps) -> float {
		 return ptr.score(bps); } )
		;
*/
    py::class_<motif_search::monte_carlo::ScorerFactory, std::shared_ptr<motif_search::monte_carlo::ScorerFactory>>(m, "MonteCarloScorerFactory")
            // ctors
            .def(py::init<>())
                    // methods
            .def("get_scorer",[] (motif_search::monte_carlo::ScorerFactory  & ptr, String const & scorer_name) -> motif_search::monte_carlo::ScorerOP {
                return ptr.get_scorer(scorer_name); } )
            ;

    py::class_<motif_search::monte_carlo::Search, std::shared_ptr<motif_search::monte_carlo::Search>>(m, "MonteCarloSearch")
            // ctors
            .def(py::init<motif_search::monte_carlo::ScorerOP, motif_search::SolutionToplogy const &,motif_search::SolutionFilterOP>())
                    // methods
            .def("clone",[] (motif_search::monte_carlo::Search const & ptr) -> motif_search::Search * {
                return ptr.clone(); } )
            .def("setup",[] (motif_search::monte_carlo::Search  & ptr, motif_search::ProblemOP p) {
                ptr.setup(p); } )
            .def("start",[] (motif_search::monte_carlo::Search  & ptr) {
                ptr.start(); } )
            .def("finished",[] (motif_search::monte_carlo::Search  & ptr) -> bool {
                return ptr.finished(); } )
            .def("next",[] (motif_search::monte_carlo::Search  & ptr) -> motif_search::SolutionOP {
                return ptr.next(); } )
            ;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_search::path_finding
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
    py::class_<motif_search::path_finding::Search, std::shared_ptr<motif_search::path_finding::Search>>(m, "PathFindingSearch")
            // ctors
            .def(py::init<motif_search::path_finding::ScorerOP,SelectorOP,motif_search::SolutionFilterOP>())
                    // methods
            .def("clone",[] (motif_search::path_finding::Search const & ptr) -> motif_search::Search * {
                return ptr.clone(); } )
            .def("setup",[] (motif_search::path_finding::Search  & ptr, motif_search:: ProblemOP prblm) {
                ptr.setup(prblm); } )
            .def("start",[] (motif_search::path_finding::Search  & ptr) {
                ptr.start(); } )
            .def("finished",[] (motif_search::path_finding::Search  & ptr) -> bool {
                return ptr.finished(); } )
            .def("next",[] (motif_search::path_finding::Search  & ptr) -> motif_search::SolutionOP {
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