#ifndef PYBIND11_THERMO_FLUCTUATION_HPP
#define PYBIND11_THERMO_FLUCTUATION_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

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


namespace thermo_fluctuation {
    namespace py = pybind11;

    void
    add_bindings(py::module_ & m ) {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// thermo_fluctuation
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // classes

        py::class_<FrameScorer, std::shared_ptr<FrameScorer>>(m, "FrameScorer")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (FrameScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
                    return ptr.score(state_1, state_2); } )
                        // inherited methods
                .def("score",[] (ThermoFlucScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
                    return ptr.score(state_1, state_2); } )
                ;

        py::class_<FrameScorerDevel, std::shared_ptr<FrameScorerDevel>>(m, "FrameScorerDevel")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (FrameScorerDevel  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
                    return ptr.score(state_1, state_2); } )
                .def("weight_d",[] (FrameScorerDevel  & ptr, float weight_d) {
                    ptr.weight_d(weight_d); } )
                .def("weight_r",[] (FrameScorerDevel  & ptr, float weight_r) {
                    ptr.weight_r(weight_r); } )
                .def("weight_d",[] (FrameScorerDevel  & ptr) ->  float {
                    return ptr.weight_d(); } )
                .def("weight_r",[] (FrameScorerDevel  & ptr) ->  float {
                    return ptr.weight_r(); } )
                        // inherited methods
                .def("score",[] (ThermoFlucScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
                    return ptr.score(state_1, state_2); } )
                ;

        py::class_<RunningAverage, std::shared_ptr<RunningAverage>>(m, "RunningAverage")
                // ctors
                .def(py::init<>())
                        // methods
                .def("Update",[] (RunningAverage  & ptr, double valIn) -> double {
                    return ptr.Update(valIn); } )
                .def("Get",[] (RunningAverage  & ptr) -> double {
                    return ptr.Get(); } )
                .def("Count",[] (RunningAverage  & ptr) -> size_t {
                    return ptr.Count(); } )
                .def("Reset",[] (RunningAverage  & ptr) {
                    ptr.Reset(); } )
                ;

        py::class_<SixDScorer, std::shared_ptr<SixDScorer>>(m, "SixDScorer")
                // ctors
                .def(py::init<String const &,structure::BasepairOP>())
        // methods
        .def("score",[] (SixDScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
            return ptr.score(state_1, state_2); } )
                // inherited methods
                .def("score",[] (ThermoFlucScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
                    return ptr.score(state_1, state_2); } )
                ;

        py::class_<ThermoFlucRelax, std::shared_ptr<ThermoFlucRelax>>(m, "ThermoFlucRelax")
                // ctors
                .def(py::init<>())
            // methods
//		.def("run_with_graph",[] (ThermoFlucRelax  & ptr) {
//		ptr.run_with_graph(); } )
                ;

        py::class_<ThermoFlucSampler, std::shared_ptr<ThermoFlucSampler>>(m, "ThermoFlucSampler")
                // ctors
                .def(py::init<>())
                        // methods
                .def("sample",[] (ThermoFlucSampler  & ptr, int steps) {
                    ptr.sample(steps); } )
                .def("setup",[] (ThermoFlucSampler  & ptr, motif_data_structure::MotifStateEnsembleTreeOP const & mset) {
                    ptr.setup(mset); } )
                .def("next",[] (ThermoFlucSampler  & ptr) -> int {
                    return ptr.next(); } )
//		.def("undo",[] (ThermoFlucSampler  & ptr) {
//		ptr.undo(); } )
                .def("to_pdb",[] (ThermoFlucSampler  & ptr, String fname, int renumber) {
                    ptr.to_pdb(fname, renumber); } )
                .def("temperature",[] (ThermoFlucSampler  & ptr) ->  float {
                    return ptr.temperature(); } )
                .def("mst",[] (ThermoFlucSampler  & ptr) ->  motif_data_structure::MotifStateTreeOP {
                    return ptr.mst(); } )
                .def("temperature",[] (ThermoFlucSampler  & ptr, float const & temp) {
                    ptr.temperature(temp); } )
                .def("randomized_start",[] (ThermoFlucSampler  & ptr, bool random_start) {
                    ptr.randomized_start(random_start); } )
                ;

        py::class_<ThermoFlucScorer, std::shared_ptr<ThermoFlucScorer>>(m, "ThermoFlucScorer")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (ThermoFlucScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
                    return ptr.score(state_1, state_2); } )
                ;

        py::class_<ThermoFlucSimulation, std::shared_ptr<ThermoFlucSimulation>>(m, "ThermoFlucSimulation")
                // ctors
                .def(py::init<>())
                        // methods
                .def("setup",[] (ThermoFlucSimulation  & ptr, motif_data_structure::MotifStateEnsembleTreeOP const & mset, int ni1, int ni2, int ei1, int ei2) {
                    ptr.setup(mset, ni1, ni2, ei1, ei2); } )
                .def("_check_sterics",[] (ThermoFlucSimulation  & ptr) ->  int {
                    return ptr._check_sterics(); } )
                .def("run",[] (ThermoFlucSimulation  & ptr) -> int {
                    return ptr.run(); } )
                .def("get_avg",[] (ThermoFlucSimulation  & ptr) -> double {
                    return ptr.get_avg(); } )
                .def("options",[] (ThermoFlucSimulation  & ptr) ->  base::Options & {
                    return ptr.options(); } )
                .def("get_int_option",[] (ThermoFlucSimulation  & ptr, String const & name) ->  float {
                    return ptr.get_int_option(name); } )
                .def("get_float_option",[] (ThermoFlucSimulation  & ptr, String const & name) ->  float {
                    return ptr.get_float_option(name); } )
                .def("get_string_option",[] (ThermoFlucSimulation  & ptr, String const & name) ->  String {
                    return ptr.get_string_option(name); } )
                .def("get_bool_option",[] (ThermoFlucSimulation  & ptr, String const & name) ->  bool {
                    return ptr.get_bool_option(name); } )
                .def("set_option_value",[] (ThermoFlucSimulation  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (ThermoFlucSimulation  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (ThermoFlucSimulation  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (ThermoFlucSimulation  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
                .def("update_var_options",[] (ThermoFlucSimulation  & ptr) {
                    ptr.update_var_options(); } )
                ;


        // classes

        py::class_<thermo_fluctuation::graph::FrameScorer, std::shared_ptr<thermo_fluctuation::graph::FrameScorer>>(m, "GraphFrameScorer")
                // ctors
                .def(py::init<>())
                        // methods
                .def("clone",[] (thermo_fluctuation::graph::FrameScorer const & ptr) -> thermo_fluctuation::graph::Scorer * {
                    return ptr.clone(); } )
                .def("score",[] (thermo_fluctuation::graph::FrameScorer  & ptr, structure::BasepairState const & state_1, structure::BasepairState const & state_2) ->  float {
                    return ptr.score(state_1, state_2); } )
                        // inherited methods
                .def("clone",[] (thermo_fluctuation::graph::Scorer const & ptr) -> thermo_fluctuation::graph::Scorer * {
                    return ptr.clone(); } )
                .def("setup",[] (thermo_fluctuation::graph::Scorer  & ptr, bool target_an_aligned_end) {
                    ptr.setup(target_an_aligned_end); } )
                .def("score",[] (thermo_fluctuation::graph::Scorer  & ptr, structure::BasepairState const & bpstate1, structure::BasepairState const & bpstate2 ) ->  float {
                    return ptr.score(bpstate1, bpstate2); } )
                ;

        py::class_<thermo_fluctuation::graph::OldFrameScorer, std::shared_ptr<thermo_fluctuation::graph::OldFrameScorer>>(m, "OldFrameScorer")
                // ctors
                .def(py::init<>())
                        // methods
                .def("clone",[] (thermo_fluctuation::graph::OldFrameScorer const & ptr) -> thermo_fluctuation::graph::Scorer * {
                    return ptr.clone(); } )
                .def("score",[] (thermo_fluctuation::graph::OldFrameScorer  & ptr, structure::BasepairState const & state_1, structure::BasepairState const & state_2) ->  float {
                    return ptr.score(state_1, state_2); } )
                        // inherited methods
                .def("clone",[] (thermo_fluctuation::graph::Scorer const & ptr) -> thermo_fluctuation::graph::Scorer * {
                    return ptr.clone(); } )
                .def("setup",[] (thermo_fluctuation::graph::Scorer  & ptr, bool target_an_aligned_end) {
                    ptr.setup(target_an_aligned_end); } )
                .def("score",[] (thermo_fluctuation::graph::Scorer  & ptr, structure::BasepairState const & bpstate1, structure::BasepairState const & bpstate2) ->  float {
                    return ptr.score(bpstate1, bpstate2); } )
                ;
/*
        py::class_<thermo_fluctuation::graph::Parameters, std::shared_ptr<thermo_fluctuation::graph::Parameters>>(m, "Parameters")
		// public attributes
		.def_readwrite("temperature", &Parameters::temperature)
		.def_readwrite("steric_radius", &Parameters::steric_radius)
		.def_readwrite("cutoff", &Parameters::cutoff)
		;
*/
        py::class_<thermo_fluctuation::graph::Sampler, std::shared_ptr<thermo_fluctuation::graph::Sampler>>(m, "Sampler")
                // ctors
                .def(py::init<motif_data_structure::MotifStateEnsembleGraph const &>())
        // methods
        .def("get_initial_state",[] (thermo_fluctuation::graph::Sampler  & ptr) -> motif_data_structure::MotifStateGraphOP {
            return ptr.get_initial_state(); } )
                .def("next",[] (thermo_fluctuation::graph::Sampler  & ptr, motif_data_structure:: MotifStateGraphOP graph) -> int {
                    return ptr.next(graph); } )
//		.def("undo",[] (Sampler  & ptr) {
//'		ptr.undo(); } )
                .def("set_temperature",[] (thermo_fluctuation::graph::Sampler  & ptr, float temp) {
                    ptr.set_temperature(temp); } )
                ;
/*
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("setup",[] (Scorer  & ptr, bool target_an_aligned_end) {
		ptr.setup(target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, structure::BasepairState const &, structure::BasepairState const &) ->  float {
		 return ptr.score(&, &); } )
		;
*/
        py::class_<thermo_fluctuation::graph::Simulation, std::shared_ptr<thermo_fluctuation::graph::Simulation>>(m, "Simulation")
                // ctors
                .def(py::init<thermo_fluctuation::graph::ScorerOP,thermo_fluctuation::graph::sterics::StericsOP>())
                        // methods
                .def("setup",[] (thermo_fluctuation::graph::Simulation  & ptr, motif_data_structure::MotifStateEnsembleGraph const & mseg, data_structure::NodeIndexandEdge const & start, data_structure::NodeIndexandEdge const & end) {
                    ptr.setup(mseg, start, end); } )
                .def("next",[] (thermo_fluctuation::graph::Simulation  & ptr) -> bool {
                    return ptr.next(); } )
                .def("write_pdbs",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) {
                    ptr.write_pdbs(name); } )
                .def("get_pdb_str",[] (thermo_fluctuation::graph::Simulation  & ptr) -> String {
                    return ptr.get_pdb_str(); } )
                .def("get_motif_graph",[] (thermo_fluctuation::graph::Simulation  & ptr) -> motif_data_structure::MotifGraphOP {
                    return ptr.get_motif_graph(); } )
                .def("get_score",[] (thermo_fluctuation::graph::Simulation  & ptr) -> float {
                    return ptr.get_score(); } )
                .def("get_int_option",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) ->  float {
                    return ptr.get_int_option(name); } )
                .def("get_float_option",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) ->  float {
                    return ptr.get_float_option(name); } )
                .def("get_string_option",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) ->  String {
                    return ptr.get_string_option(name); } )
                .def("get_bool_option",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) ->  bool {
                    return ptr.get_bool_option(name); } )
                .def("set_option_value",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
                .def("set_option_value",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
                ;

        py::class_<thermo_fluctuation::graph::logging::TargetBPInfoLogger, std::shared_ptr<thermo_fluctuation::graph::logging::TargetBPInfoLogger>>(m, "TargetBPInfoLogger")
                // ctors
                .def(py::init<String const &>())
        // methods
        .def("clone",[] (thermo_fluctuation::graph::logging::TargetBPInfoLogger const & ptr) -> thermo_fluctuation::graph::logging::Logger * {
            return ptr.clone(); } )
                .def("setup",[] (thermo_fluctuation::graph::logging::TargetBPInfoLogger  & ptr, motif_data_structure::MotifStateGraphOP msg, data_structure::NodeIndexandEdge const & start, data_structure::NodeIndexandEdge const & end) {
                    ptr.setup(msg, start, end); } )
                .def("log",[] (thermo_fluctuation::graph::logging::TargetBPInfoLogger  & ptr, motif_data_structure::MotifStateGraphOP msg, float score) {
                    ptr.log(msg, score); } )
                        // inherited methods
                .def("clone",[] (thermo_fluctuation::graph::logging::Logger const & ptr) -> thermo_fluctuation::graph::logging::Logger * {
                    return ptr.clone(); } )
                .def("setup",[] (thermo_fluctuation::graph::logging::Logger  & ptr, motif_data_structure:: MotifStateGraphOP graph, data_structure::NodeIndexandEdge const & nandi1, data_structure::NodeIndexandEdge const &nandi2) {
                    ptr.setup(graph, nandi1, nandi2); } )
                .def("log",[] (thermo_fluctuation::graph::logging::Logger  & ptr, motif_data_structure:: MotifStateGraphOP graph, float  value) {
                    ptr.log(graph, value); } )
                ;


        // classes

        py::class_<thermo_fluctuation::graph::sterics::NoSterics, std::shared_ptr<thermo_fluctuation::graph::sterics::NoSterics>>(m, "NoSterics")
                // ctors
                .def(py::init<>())
                        // methods
                .def("clone",[] (thermo_fluctuation::graph::sterics::NoSterics const & ptr) -> thermo_fluctuation::graph::sterics::Sterics * {
                    return ptr.clone(); } )
                .def("clash",[] (thermo_fluctuation::graph::sterics::NoSterics  & ptr, motif_data_structure::MotifStateGraphOP msg) -> bool {
                    return ptr.clash(msg); } )
                        // inherited methods
                .def("clone",[] (thermo_fluctuation::graph::sterics::Sterics const & ptr) -> thermo_fluctuation::graph::sterics::Sterics * {
                    return ptr.clone(); } )
                .def("clash",[] (thermo_fluctuation::graph::sterics::Sterics  & ptr, motif_data_structure:: MotifStateGraphOP graph) -> bool {
                    return ptr.clash(graph); } )
                ;

        py::class_<thermo_fluctuation::graph::sterics::SelectiveSterics, std::shared_ptr<thermo_fluctuation::graph::sterics::SelectiveSterics>>(m, "SelectiveSterics")
                // ctors
                .def(py::init<Ints const,Ints const,float>())
        // methods
        .def("clone",[] (thermo_fluctuation::graph::sterics::SelectiveSterics const & ptr) -> thermo_fluctuation::graph::sterics::Sterics * {
            return ptr.clone(); } )
                .def("clash",[] (thermo_fluctuation::graph::sterics::SelectiveSterics  & ptr, motif_data_structure::MotifStateGraphOP msg) -> bool {
                    return ptr.clash(msg); } )
                        // inherited methods
                .def("clone",[] (thermo_fluctuation::graph::sterics::Sterics const & ptr) -> thermo_fluctuation::graph::sterics::Sterics * {
                    return ptr.clone(); } )
                .def("clash",[] (thermo_fluctuation::graph::sterics::Sterics  & ptr, motif_data_structure:: MotifStateGraphOP graph) -> bool {
                    return ptr.clash(graph); } )
                ;
/*
        py::class_<Sterics, std::shared_ptr<Sterics>>(m, "Sterics")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Sterics const & ptr) -> Sterics * {
		 return ptr.clone(); } )
		.def("clash",[] (Sterics  & ptr, motif_data_structure:: MotifStateGraphOP graph) -> bool {
		 return ptr.clash(graph); } )
		;
*/
    }
}

#endif // PYBIND11_THERMO_FLUCTUATION_HPP