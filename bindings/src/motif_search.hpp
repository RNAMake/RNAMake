#ifndef PYBIND11_MOTIF_SEARCH_HPP
#define PYBIND11_MOTIF_SEARCH_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>
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


namespace motif_search::monte_carlo {
    namespace py = pybind11;
    void
    add_bindings(py::module_ & );
}

namespace motif_search::exhaustive {
    namespace py = pybind11;
    void
    add_bindings(py::module_ & );
}

namespace motif_search::path_finding {
    namespace py = pybind11;
    void
    add_bindings(py::module_ & );
}

namespace motif_search {
    namespace py = pybind11;

    void
    add_bindings(py::module_ & m ) {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// motif_search
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        auto monte_carlo = m.def_submodule("monte_carlo");
        monte_carlo::add_bindings(monte_carlo);

        auto exhaustive = m.def_submodule("exhaustive");
        exhaustive::add_bindings(exhaustive);

        auto path_finding = m.def_submodule("path_finding");
        path_finding::add_bindings(path_finding);

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


    }

}

namespace motif_search::monte_carlo {
    namespace py = pybind11;
    void
    add_bindings(py::module_ &  m) {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// motif_search::monte_carlo
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // classes

        py::class_<motif_search::monte_carlo::DefaultScorer, std::shared_ptr<motif_search::monte_carlo::DefaultScorer>>(m, "DefaultScorer")
                // ctors
                .def(py::init<>())
                        // methods
                .def("clone",[] (motif_search::monte_carlo::DefaultScorer const & ptr) -> motif_search::monte_carlo::Scorer * {
                    return ptr.clone(); } )
                .def("score",[] (motif_search::monte_carlo::DefaultScorer  & ptr, structure::BasepairState const & bps) ->  float {
                    return ptr.score(bps); } )
                        // inherited methods
                .def("clone",[] (motif_search::monte_carlo::Scorer const & ptr) -> motif_search::monte_carlo::Scorer * {
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
        py::class_<motif_search::monte_carlo::ScorerFactory, std::shared_ptr<motif_search::monte_carlo::ScorerFactory>>(m, "Factory")
                // ctors
                .def(py::init<>())
                        // methods
                .def("get_scorer",[] (motif_search::monte_carlo::ScorerFactory  & ptr, String const & scorer_name) -> motif_search::monte_carlo::ScorerOP {
                    return ptr.get_scorer(scorer_name); } )
                ;

        py::class_<motif_search::monte_carlo::Search, std::shared_ptr<motif_search::monte_carlo::Search>>(m, "Search")
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
    }
}

namespace motif_search::exhaustive {
    namespace py = pybind11;
    void
    add_bindings(py::module_ & m ) {

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// motif_search::exhaustive
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // classes

        py::class_<motif_search::exhaustive::DefaultScorer, std::shared_ptr<motif_search::exhaustive::DefaultScorer>>(m, "DefaultScorer")
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
        py::class_<motif_search::exhaustive::ScorerFactory, std::shared_ptr<motif_search::exhaustive::ScorerFactory>>(m, "ScorerFactory")
                // ctors
                .def(py::init<>())
                        // methods
                .def("get_scorer",[] (motif_search::exhaustive::ScorerFactory  & ptr, String const & scorer_name) -> motif_search::exhaustive::ScorerOP {
                    return ptr.get_scorer(scorer_name); } )
                ;

        py::class_<motif_search::exhaustive::Search, std::shared_ptr<motif_search::exhaustive::Search>>(m, "Search")
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



    }
}

namespace motif_search::path_finding {
    namespace py = pybind11;
    void
    add_bindings(py::module_ & m ) {
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
                .def("clone",[] (AstarScorer const & ptr) -> motif_search::path_finding::Scorer * {
                    return ptr.clone(); } )
                .def("score",[] (AstarScorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.score(node); } )
                .def("score",[] (AstarScorer  & ptr, motif::MotifState & ms, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.score(ms, node); } )
                        // inherited methods
                .def("clone",[] (motif_search::path_finding::Scorer const & ptr) -> motif_search::path_finding::Scorer * {
                    return ptr.clone(); } )
                .def("set_target",[] (motif_search::path_finding::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                    ptr.set_target(target, target_an_aligned_end); } )
                .def("score",[] (motif_search::path_finding::Scorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.score(node); } )
                .def("score",[] (motif_search::path_finding::Scorer  & ptr, motif::MotifState & state, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.score(state, node); } )
                .def("accept_score",[] (motif_search::path_finding::Scorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.accept_score(node); } )
                .def("set_dummy",[] (motif_search::path_finding::Scorer  & ptr, float dummy) {
                    ptr.set_dummy(dummy); } )
                ;

        py::class_<GreedyScorer, std::shared_ptr<GreedyScorer>>(m, "GreedyScorer")
                // ctors
                .def(py::init<>())
                        // methods
                .def("clone",[] (GreedyScorer const & ptr) -> motif_search::path_finding::Scorer * {
                    return ptr.clone(); } )
                .def("score",[] (GreedyScorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.score(node); } )
                .def("score",[] (GreedyScorer  & ptr, motif::MotifState & ms, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.score(ms, node); } )
                        // inherited methods
                .def("clone",[] (motif_search::path_finding::Scorer const & ptr) -> motif_search::path_finding::Scorer * {
                    return ptr.clone(); } )
                .def("set_target",[] (motif_search::path_finding::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                    ptr.set_target(target, target_an_aligned_end); } )
                .def("score",[] (motif_search::path_finding::Scorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.score(node); } )
                .def("score",[] (motif_search::path_finding::Scorer  & ptr, motif::MotifState & state, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.score(state, node); } )
                .def("accept_score",[] (motif_search::path_finding::Scorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                    return ptr.accept_score(node); } )
                .def("set_dummy",[] (motif_search::path_finding::Scorer  & ptr, float dummy) {
                    ptr.set_dummy(dummy); } )
                ;

        py::class_<motif_search::path_finding::Node, std::shared_ptr<motif_search::path_finding::Node>>(m, "Node")
                // ctors
                .def(py::init<motif::MotifStateOP,motif_search::path_finding::NodeOP,float,int,int,int>())
                        // methods
                .def("level",[] (motif_search::path_finding::Node const & ptr) -> int {
                    return ptr.level(); } )
                .def("size",[] (motif_search::path_finding::Node const & ptr) ->  int {
                    return ptr.size(); } )
                .def("node_type",[] (motif_search::path_finding::Node const & ptr) ->  int {
                    return ptr.node_type(); } )
                .def("ss_score",[] (motif_search::path_finding::Node const & ptr) ->  float {
                    return ptr.ss_score(); } )
                .def("score",[] (motif_search::path_finding::Node const & ptr) ->  float {
                    return ptr.score(); } )
                .def("parent_end_index",[] (motif_search::path_finding::Node const & ptr) ->  int {
                    return ptr.parent_end_index(); } )
                .def("state",[] (motif_search::path_finding::Node const & ptr) ->  motif::MotifStateOP {
                    return ptr.state(); } )
                .def("parent",[] (motif_search::path_finding::Node const & ptr) ->  motif_search::path_finding::NodeOP {
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
        py::class_<motif_search::path_finding::Search, std::shared_ptr<motif_search::path_finding::Search>>(m, "Search")
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
}
#endif // PYBIND11_MOTIF_SEARCH_HPP